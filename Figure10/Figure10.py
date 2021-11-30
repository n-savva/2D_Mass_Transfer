# -*- coding: utf-8 -*-
"""
Data generation for Figure 10
"""
import sys
sys.path.append("../main/")
from ODEdrop2D import *
from pdeloader import *
from matplotlib.gridspec import GridSpec
from scipy.io import loadmat

V = lambda t: Vperiodic(t,m=20,p=P,Vavg=2.5,Vamp=1.5)
theta = lambda x: 1 + 0.1*np.cos(1.6*np.pi*x) + 0.3*np.sin(10.*np.pi*x/3.)


fig = plt.figure(figsize=(8,2.5))
gs = GridSpec(1, 2, figure=fig,wspace=0.24)

for i,P,T,label in zip([0,1],[200,400],[800,2000],['a','b']):
    ax = fig.add_subplot(gs[i],adjustable='box')
    
    # Plot the PDE data
    t_pde, __ , l_pde = pdeloader('Figure10_PDE_{0}.mat'.format(P))
    ax.plot(t_pde[::20],l_pde[::20],'k',lw=0.25,label='PDE')

    # Solve the ODE problem
    drop = ODEdrop2D(ic=(1,-1),t_end=T,het=theta,V=V)
    t = np.linspace(0,T,400)
    drop.solve()
    ax.plot(t,np.mean(drop.evaluate(t),axis=0),'--',lw=1,label='ODEs')
    ax.set_xlim((0,T))
    ax.set_xlabel('$t$')
    ax.set_ylabel('$\\ell$')
    ax.set_title('$p={0}$'.format(P))
    ax.text(-0.22*T,ax.get_ylim()[1],'({0})'.format(label))

ax.legend()
plt.savefig('Figure10.png', bbox_inches='tight',dpi=200)