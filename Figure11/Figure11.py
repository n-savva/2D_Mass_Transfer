# -*- coding: utf-8 -*-
"""
Data generation for Figure 11
"""
import sys
sys.path.append("../main/")
from ODEdrop2D import *
from pdeloader import *
from matplotlib.gridspec import GridSpec
from scipy.io import loadmat
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle
from scipy.interpolate import interp1d


def V(t):
    return 2 - w*t, -w

data = loadmat('Figure11_Het.mat')
g = interp1d(data['x'].flatten(),data['G'].flatten())

# Array of w's
W = np.array([0.025,0.005,0.001])
Tf = 2/W-1 

Ts = (np.hstack((np.arange(50,350,50),[395])),np.hstack((np.arange(250,1500,250),[1990])))

fig = plt.figure(figsize=(8,6))
gs = GridSpec(2, 2, figure=fig,wspace=0.24) 

for iplot,w,t_plot,label in zip((0,1),(0.005,0.001),Ts,('a','b')):
    ax = fig.add_subplot(gs[0,iplot],adjustable='box')
    plt.sca(ax)
    PDEdata = loadmat("Figure11_PDE{0}.mat".format(int(1e3*w)))
    t_pde = PDEdata['tPDE'].flatten()

    drop = ODEdrop2D(ic=(1,-1),t_end=2/w-1,het=g,V=V,slip=1e-4)
    drop.solve()
    
    for ti in t_plot:
        i = np.argmin(np.abs(t_pde-ti))
        h = PDEdata['h'][i,:]
        x = PDEdata['xp'][i,:]
        plt.plot(x,h,'k',lw=0.25)
                   
    drop.draw(t_plot,ls='--',color='tab:blue',lw=1,xlim=(-2.5,2.5))
    ax.text(-3.5,ax.get_ylim()[1],'({0})'.format(label))
    ax.set_xlabel('$x$')
    ax.set_ylabel('$h$')
    
ax = fig.add_subplot(gs[1,0],adjustable='box')
ax.plot(data['x'].flatten(),data['G'].flatten(),'k',lw=0.5)
ax.plot(data['x'].flatten(),1-0.2*np.tanh(50*np.cos(np.pi*data['x'].flatten())),'--',lw=1,color='tab:red')
ax.set_xlim((-2.5,2.5))
ax.set_xlabel('$x$')
ax.set_ylabel('$\\theta$')
ax.text(-3.5,ax.get_ylim()[1],'(c)')

ax = fig.add_subplot(gs[1,1],adjustable='box')
ax.set_xlim((0,1))
plt.sca(ax)

for i in range(3):
    w = W[i]
    
    PDEdata = loadmat("Figure11_PDE{0}.mat".format(int(1e3*w)))
    l_pde = 0.5*(PDEdata['a'].flatten()+PDEdata['b'].flatten())
    t_pde = PDEdata['tPDE'].flatten()

    # Solve the ODE problem
    drop = ODEdrop2D(ic=(1,-1),t_end=Tf[i],het=g,V=V,flux=None,slip=1e-4)
    drop.solve()
    
    t = np.linspace(0,Tf[i],500)
    ab  = drop.evaluate(t)
    l = 0.5*(ab[0]+ab[1])
    plt.plot(t_pde[::5]/t_pde[-1],l_pde[::5],'k',lw=0.25)
    plt.plot(t/Tf[i],l,'--',lw=1,label='$w = {0}$'.format(w))
    
ax.text(-0.2,ax.get_ylim()[1],'(d)')
ax.set_ylabel('$\ell$')
ax.set_xlabel('$t/T_f$')
plt.legend()
plt.savefig('Figure11.png', bbox_inches='tight',dpi=200)