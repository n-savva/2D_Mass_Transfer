# -*- coding: utf-8 -*-
"""
Data generation for Figure 09
"""
import sys
sys.path.append("../main/")
from ODEdrop2D import *
from pdeloader import *
from matplotlib.gridspec import GridSpec

V = lambda t: Vperiodic(t,m=20,p=100,Vavg=2,Vamp=1)
theta = lambda x: 1 + .1*np.cos(1.6*np.pi*x) + 0.2*np.sin(0.2*np.pi*x)

# Solve the ODE problem
drop = ODEdrop2D(ic=(1,-1),t_end=1000,het=theta,V=V)
t = np.arange(0,1001)

fig = plt.figure(figsize=(8,3))
gs = GridSpec(1, 2, figure=fig,wspace=0.24)
ax0 = fig.add_subplot(gs[0],adjustable='box')
ax1 = fig.add_subplot(gs[1],adjustable='box')

# Load PDE
t_pde,d_pde , l_pde = pdeloader('Figure09_PDE.mat')

# Solve ODE
drop.solve()
ab = drop.evaluate(t)

# Plot of the midpoint
ax0.plot(V(t_pde)[0],l_pde,'k',lw=0.25)
ax0.plot(V(t)[0],np.mean(ab,axis=0),'--',lw=1)
ax0.set_xlabel('$v$')
ax0.set_ylabel('$\\ell$')
ax0.text(0.3,ax0.get_ylim()[1],'(a)')

# Plot of the radius
ax1.plot(V(t_pde)[0],d_pde,'k',lw=0.25,label='PDE')
ax1.plot(V(t)[0],0.5*(ab[0]-ab[1]),'--',lw=1,label='ODEs')
ax1.set_xlabel('$v$')
ax1.set_ylabel('$d$')
ax1.legend()
ax1.text(0.5,ax1.get_ylim()[1],'(b)')

plt.savefig('Figure09.png', bbox_inches='tight',dpi=200)