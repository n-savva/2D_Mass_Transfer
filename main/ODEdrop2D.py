# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from tqdm import tqdm
from matplotlib.colors import LinearSegmentedColormap

# Periodic volume variations
def Vperiodic(t,m=20,p=100,Vavg=2,Vamp=1):
    a = 2*np.pi/p
    A = np.sqrt(1+(m*np.cos(a*t))**2)
    V = Vavg + Vamp*np.arctan(m*np.sin(a*t)/A)/np.arctan(m)
    Vdot = Vamp*a*m*np.cos(a*t)/A/np.arctan(m)
    return V,Vdot

class ODEdrop2D(object):
    """
    A class for a drop object.
    
    Attributes
    ----------
    slip: float  [1e-4]
        Slip length.
    V: float [2]
        Droplet area. It can be a constant or a function of time. If it is a
        function of time, the function should also return its tme derivative.
    n: int [100]
        Number of quadrature points for the integral terms.
    het: float or callable [1.0]
        The heterogeneity profile.
    ic: tuple [(1.0,-1.0)]
        The initial condition for the right and left contact points
    t_end: float [100]
        The final time.
    flux: Either None or a callable [None]
        None: appropriate for parabolic flux or for constant area
        If function of time only, it corresponds to the delta-localized flux.
        If function of x and t, allows for localized fluxes. Output is
        normalized so that its integral is equal to unity.
    method: str ['RK45']
        The method to be used with solve_ivp. You may also use 'LSODA'.
    soltion: OdeSolution
        OdeSolution instance containing the solution to the equations; See 
        documentation for solve_ivp.
    
    Methods
    -------
    ic(ic), het(het), V(vol) 
        Sets the attributes for ic het and V
    
    solve()
        Computes the solution if t_end is specified.
        
    drawcl(T,color='b',style='-')    
        Draws contact line shapes for given values of time in T.

        Parameters
        ----------
        T : array_like
            The times at which the shapes are to be plotted.
        color: string, optional
            The color of the plot. The default is 'b'.
        style: string, optional
            The line style for the plot. The default is '-'
        
        Raises
        ------
        Exception
            - If t_end is undefined.
            - If solution has not yet been computed.
            - If some element of T lies outside the time range.

        Returns
        -------
        None

    
    angle(t)    
        Returns an array with the apparent angle at a specified time, t.

        Parameters
        ----------
        t : float 
           The time for which the angle is to be returned. The default is None.

        Raises
        ------
        Exception
            - If t lies outside the solution range
            - If solution has yet to be computed
            - If t is not a scalar

        Returns
        -------
        angle: array_like
            The apparent contact angle. 

        
    evaluate(t)
        Returns the contact point positions at prescribed times.

        Parameters
        ----------
        t : float or array_like
            The times at which the contact points are to be returned. 
            An exceptionis thrown if some element of t lies outside the solution range.
            
        Returns
        -------
        a,b: array_like
            Right and left contact points
        
    
    draw(t,kind='midpoint',*args,**kwargs)    
        Makes plots of the generated data

        Parameters
        ----------
        t : float or array_like
            Time array for plotting. The default is None.
        kind : the kind of data to be plotted, optional
            It can be either 'midpoint', 'width', 'left' or 'right'. The 
            default is 'midpoint'.
        **args,**kwargs : Additional (optional) arguments to be passed on to
            the plotting routine.
    
    """
    def __init__(self, slip = 1e-4, V=2, n = 100, het = 1., 
                 ic = (1.0,-1.0), flux=None, t_end=100,method='RK45'):
        
        # Discretization
        self.slip = slip
        self.n = n
        self.__p = np.log(0.5*slip)+1.0
        self.__s = np.array([1,-1])

        # Parse simulation parameters
        self.V = V
        self.flux = flux
        
        if flux is not None:
            self.__isdelta = True
            if callable(flux):
                self.__isfunc = True
                if flux.__code__.co_argcount==2:
                    self.__isdelta = False
                            
        # Initial condition
        self.ic = ic
        
        # ODE integrator
        self.t_end = t_end
        self.method = method
        self.solution = None
        
        # Chemical Heterogeneity
        self.het = het
        self.X, self.W = np.polynomial.legendre.leggauss(self.n)
        self.logx = self.W*(0.5*np.log((1+self.X)/(1-self.X)) \
            + self.__s[:,None]/(1-self.__s[:,None]*self.X))
                         
    # Volume property
    @property
    def V(self):
        return self._V
    
    @V.setter
    def V(self,value):
        if not callable(value):
            self._V = lambda t: (value,0)
        else:
            self._V = value 
                                
    # Intial Condition
    @property
    def ic(self):
        return self._ic
    
    @ic.setter
    def ic(self,value):
        self._ic = value
        self.solution = None
   
    # Heterogeneity profile
    @property
    def het(self):
        return self._g
    
    @het.setter
    def het(self,value):
        if not callable(value):
            self._g = lambda x: np.full(x.shape,value,dtype='float64')
        else:
            self._g = value
        self.solution = None
   
    # Evaluate Radius via BIM
    def __angle(self,Vo,Ro): 
        return 1.5*Vo/Ro**2

    # Progress Bar Display    
    def __pbar(self,t,y,pbar,state):
        last_t, dt = state
        n = int((t-last_t)/dt)
        pbar.update(n)
        pbar.set_description("t = %1.2e" % t)
        state[0] = last_t + dt*n
        return 1   

    # ODE
    def __ode(self,t,U,pbar,state):
        # Centroid and harmonics
        d = 0.5*(U[0]-U[1])
        V,Vdot = self.V(t)
        
        # Local and apparent contact angles
        θs = self._g(U)
        θ = self.__angle(V, d)

        # Flux term
        Flux = 0
        if self.flux is not None:
            if self.__isdelta:
                if self.__isfunc:
                    xo = self.flux(t)
                else:
                    xo = self.flux
                 
                I = 0.5/d*np.log((xo - U[1])/(U[0] - xo)) + 1./(U-xo)
            else:
                q = self.flux(t,d*self.X + (U[0]+U[1])*0.5)       
            
                I = np.dot(self.logx,q)/np.dot(self.W,q)/d
                
            Flux = self.__s * I*Vdot/θ - 1.5*Vdot/d/θ
                
        
        K = (θ**3 - θs**3)/3. + Flux
        L = np.log(d*θs) - self.__p
        
        return self.__s*(K*L[::-1] + K[::-1])/(L[0]*L[1]-1.)
        
        
    # Solve Method
    def solve(self):
        """
        Solves the system.

        """
        if self.t_end is None:
            raise Exception("Undefined t_end")
        else:
           
            print("\nSolving until t = %1.2f\n" % (self.t_end),end='',flush=True)
            with tqdm(total=100,unit="%") as pbar:
                self.solution = solve_ivp(self.__ode, (0,self.t_end), 
                                          np.array(self.ic), 
                                          method=self.method,
                                          dense_output=True,
                                          events=self.__pbar,
                                          atol=1e-8,rtol=1e-8,
                                          args=[pbar,[0,self.t_end/100]])
                
    
    # Method for returning the contact point positions
    def evaluate(self,t=None):
        """
        Returns the contact point positions at prescribed times.

        Parameters
        ----------
        t : float or array_like
            The times at which the contact points are to be returned. 
            An exceptionis thrown if some element of t lies outside the 
            solution range.
            
        Returns
        -------
        a,b: array_like
            Right and left contact points 
            
        """
        if t is not None and self.solution is not None:
            if not isinstance(t, (list, tuple, np.ndarray)):
                t = [t]
            L = len(t)
            a = np.zeros(L)
            b = np.zeros(L)
            
            for i in range(L):
                if t[i]<0 or t[i]>self.t_end:
                    raise Exception('Time out of range')
                U = self.solution.sol(t[i])
                a[i],b[i] = U[0],U[1]
                     
            if L==1:
                a, b = a[0], b[0]
    
            return a, b
    
    def plot(self,t=None,kind='midpoint',**kwargs):
        """
        Makes plots of the generated data

        Parameters
        ----------
        t : float or array_like
            Time array for plotting. The default is None.
        kind : the kind of data to be plotted, optional
            It can be either 'midpoint', 'width', 'left', 'right' or 'angle'.
            The default is 'midpoint'.
        **kwargs : Additional (optional) arguments to be passed on to
            the plotting routine.

        """
        if self.solution is not None and t is not None:
            
            if kind.lower() == 'angle':
                data = self.angle(t)
            else:
                ab = self.evaluate(t)
            
                if kind.lower()=='midpoint':
                    data = 0.5*(ab[0] + ab[1])
                elif kind.lower() =='width':
                    data = 0.5*(ab[0] - ab[1])
                elif kind.lower()=='right':
                    data = ab[0]
                elif kind.lower()=='left':
                    data = ab[1]
                else:
                    raise Exception('Unknown kind of data selected')
               
            plt.plot(t,data,**kwargs)
            
        
    # Draw droplet
    def draw(self,t=None,xlim=(-1.5,1.5),**kwargs):
        """
        Draws droplet and surface profile at prescribed times.

        Parameters
        ----------
        t : float or array_like
            Time array for plotting. The default is None.
        xlim : tuple of floats, optional
            The limit of the x-axis for drawing the surface. The default is 
            (-1.5,1.5).
        **kwargs  : Additional (optional) arguments to be passed on to
            the plotting routine.

        """
        if self.solution is not None and t is not None:
            if not isinstance(t, (list, tuple, np.ndarray)):
                t = [t]
            
            x = np.linspace(-1,1)
            for ti in t:
                ab = self.evaluate(ti)
                
                plt.plot(0.5*((ab[0]-ab[1])*x+ab[0]+ab[1]),\
                         1.5*self.V(ti)[0]*(1-x**2)/(ab[0]-ab[1]),**kwargs)
    
            x = np.linspace(xlim[0],xlim[1],num=300)
            plt.plot(xlim,(0,0),'k',lw=2)
            plt.pcolor(x,[-.05,0],np.array([self._g(x)])[[0],:-1],\
                       cmap=LinearSegmentedColormap.from_list('custom greys', [(1,1,1),(0.5,0.5,0.5)], N=256),\
                       shading='flat')
                
    # Method for returning the contact point positions
    def angle(self,t=None):
        """
        Returns the apparent angle for prescribed times.

        Parameters
        ----------
        t : float or array_like
            The times at which the contact points are to be returned. 
            An exception is thrown if some element of t lies outside the solution range.
            
        Returns
        -------
        theta: array_like
            Apparent contact angle
        """
        if t is not None and self.solution is not None:
            if not isinstance(t, (list, tuple, np.ndarray)):
                t = [t]
            L = len(t)
            theta = np.zeros(L)
            
            for i in range(L):
                if t[i]<0 or t[i]>self.t_end:
                    raise Exception('Time out of range')
                U = self.solution.sol(t[i])
                V,__ = self.V(t[i])
                theta[i] = self.__angle(V,0.5*(U[0]-U[1]))
                     
            if L==1:
                theta = theta[0]
    
            return theta