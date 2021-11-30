# -*- coding: utf-8 -*-
"""
Data generation for Figure 07
"""
import sys
sys.path.append("../main/")
from ODEdrop2D import *
from pdeloader import *
from matplotlib.gridspec import GridSpec
from scipy.io import loadmat
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle

def q(t,x):
    S = 20
    return np.exp(-S*(x-xo)**2)
    

def V(t):
    return  Vperiodic(t,m=20,p=250,Vavg=3,Vamp=1.25)

g = lambda x: 1 + .1*np.tanh(5*np.cos(np.pi*x))

# Solve the ODE problem
drop = ODEdrop2D(ic=(1,-1),t_end=1000,het=g,V=V,flux=q)
t = np.linspace(0,1000,201)

fig = plt.figure(figsize=(8,2.5))
gs = GridSpec(1, 2, figure=fig,wspace=0.24)
ax0 = fig.add_subplot(gs[0],adjustable='box')
ax1 = fig.add_subplot(gs[1],adjustable='box')


for xo in [0.25,0.75]:
    t_pde, d_pde, l_pde = pdeloader("Figure07_PDE_0{0}.mat".format(int(100*xo)))
    
    drop.solve()
    ab  = drop.evaluate(t)
    d = 0.5*(ab[0]-ab[1])
    l = 0.5*(ab[0]+ab[1])

    
    ax0.plot(t_pde[::5],l_pde[::5],'k',lw=0.25)
    ax1.plot(t_pde[::5],d_pde[::5],'k',lw=0.25)
    ax1.plot(t,d,'--',label='$x_0={0}$'.format(xo))
    ax0.plot(t,l,'--')
   
ax1.legend()

for ax,ylab,label in zip([ax0,ax1],['$\\ell$','$d$'],['a','b']):
    ax.set_xlim((0,1000))
    ax.set_xlabel('$t$')
    ax.set_ylabel(ylab)
    ax.text(-200,ax.get_ylim()[1],'({0})'.format(label))

plt.savefig('Figure07.png', bbox_inches='tight',dpi=200)