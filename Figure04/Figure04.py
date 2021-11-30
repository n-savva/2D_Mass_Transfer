# -*- coding: utf-8 -*-
"""
Data generation for Figure 04
"""
import sys
sys.path.append("../main/")
from ODEdrop2D import *
from pdeloader import *
from matplotlib.gridspec import GridSpec
from scipy.io import loadmat

def V(t):
    return  Vperiodic(t,m=20,p=200,Vavg=1,Vamp=0.5)

g = lambda x: 1 + .5*np.tanh(30*np.cos(c*x))

fig = plt.figure(figsize=(6,3),tight_layout=True)
gs = GridSpec(2, 2, figure=fig)


for i,c,label,lims in zip(range(4),[9,10,12,13],['a','b','c','d'],\
                          [(1.05,1.25),(1.05,1.25),(1.30,1.46),(1.30,1.46)]):
    # Create axes
    ax = fig.add_subplot(gs[i],adjustable='box')
    
    # Load PDE data
    t_pde, d_pde, l_pde = pdeloader("Figure04_PDE{0}.mat".format(label))
    
    # Solve ODE system
    drop = ODEdrop2D(ic=(1.,-1),t_end=600,het=g,V=V)
    drop.solve()
    sol = drop.evaluate(t_pde[::20])
    
    # Plots
    plt.sca(ax)
    ax.plot(t_pde[::20],d_pde[::20],'k',lw=0.25,label='PDE')
    drop.plot(t_pde[::20],kind='width',ls='--',lw=1,label='ODEs')
    ax.set_xlim((200,400))
    ax.set_ylim(lims)
    
    # labels
    if i%2==0:
        ax.set_ylabel('$d$')
        
    if i>1:
        ax.set_xlabel('$t$')
    
    if i==3:
        plt.legend()
    
    ax.text(150,ax.get_ylim()[1],'({0})'.format(label))

plt.savefig('Figure04.png',bbox_inches='tight', dpi=200)