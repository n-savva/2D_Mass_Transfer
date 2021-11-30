# -*- coding: utf-8 -*-
"""
Data generation for Figure 08
"""
import sys
sys.path.append("../main/")
from ODEdrop2D import *
from pdeloader import *
from matplotlib.gridspec import GridSpec
from scipy.io import loadmat
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle

def q(t):
    n = 1 + (t//100) # Cycle
    
    if t < (100*n-50):
        return 0.6 + 0.4*n
    else:
        return -1 + 0.4*n

V = lambda t: Vperiodic(t-25,m=20,p=100,Vavg=2,Vamp=1)


# Solve the ODE problem
drop = ODEdrop2D(ic=(np.sqrt(1.5),-np.sqrt(1.5)),t_end=600,V=V,flux=q)
drop.solve()
t = np.linspace(0,600,301)
t_in = sum([[0+ti,ti+49.999,np.nan] for ti in np.arange(0,600,100)],[])
t_out = sum([[50+ti,ti+99.999,np.nan] for ti in np.arange(0,600,100)],[])


ab = drop.evaluate(t)


fig = plt.figure(figsize=(8,4))
gs = GridSpec(1, 2, figure=fig,wspace=0.25)

ax0 = fig.add_subplot(gs[0],adjustable='box')

T,X = np.meshgrid(t,[-1,0,1])
U = np.tile(0.5-0.2*np.sign(V(t)[1]),(3,1))

ax0.pcolor(0.5*(ab[0]-ab[1])*X+0.5*(ab[0]+ab[1]),T,U[:-1,:-1],\
           cmap='Greys',vmin=0,vmax=1,shading='flat')

ax0.plot(ab[0],t,ab[1],t)
ax0.plot([q(ti) for ti in t_in],t_in,'w',lw=2)
ax0.plot([q(ti) for ti in t_out],t_out,'w--',lw=2)


ax1 = fig.add_subplot(gs[1],adjustable='box')


PDEdata = loadmat('Figure08_PDEb.mat')
a_pde = PDEdata['a'].flatten()
b_pde = PDEdata['b'].flatten()
t_pde = PDEdata['t_final_PDE'].flatten()

ax1.plot(t_pde,0.5*(a_pde+b_pde),'k',lw=0.25)
ax1.plot(t,0.5*(ab[0]+ab[1]),'--',lw=1)

ax1.text(-150,ax1.get_ylim()[1],'(b)')
ax1.set_ylabel('$\\ell$')
ax1.set_xlabel('$t$')

ax0.text(-3.25,ax0.get_ylim()[1],'(a)')
ax0.set_ylabel('$t$')
ax0.set_xlabel('$x$')

plt.savefig('Figure08.png', bbox_inches='tight',dpi=200)