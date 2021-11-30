from scipy.io import loadmat
import numpy as np

def pdeloader(file,t=None,V = None):
    PDEdata = loadmat(file)
    a_left = PDEdata['yPDE'][:,-1]
    a_right = PDEdata['yPDE'][:,-2]
    d_pde = 0.5*(a_right - a_left)
    l_pde = 0.5*(a_right + a_left)
    t_pde = PDEdata['tPDE'].flatten()
    
    if t is None or V is None:
        return t_pde, d_pde, l_pde
    else:
        D1y = PDEdata['D1y']       # Differentiation matrix
        y_pde = PDEdata['yPDE']    # PDE solution
        y = PDEdata['y'].flatten() # Grid
        
        L = []        
        if not isinstance(t, (list, tuple, np.ndarray)):
            t = [t]
            
        for ti in t:
            i = np.argmin(np.abs(t_pde-ti))
            Hin =  D1y@y_pde[i,:-2]  + 0.75*V(t_pde[i])[0]*(1-y[1:-1]**2)/d_pde[i]
            L.append((d_pde[i]*y+l_pde[i],np.concatenate(([0],Hin,[0]))))
                
        return (t_pde,d_pde,l_pde), L
        
        
        