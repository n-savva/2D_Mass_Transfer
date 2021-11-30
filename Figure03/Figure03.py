# -*- coding: utf-8 -*-
"""
Data generation for figure 03
"""
import sys
sys.path.append("../main/")
from ODEdrop2D import *
from pdeloader import *
from matplotlib.gridspec import GridSpec

def V(t):
    return  Vperiodic(t,m=20,p=200,Vavg=1,Vamp=0.5)

g = lambda x: 1 + .5*np.tanh(30*np.cos(11*x))

# LOAD PDE DATA
t_plots = np.arange(50,150,10)
(t_pde,d_pde,l_pde), L = pdeloader("Figure03_PDE.mat",t=t_plots,V=V)
theta_pde = 1.5*V(t_pde)[0]/d_pde**2

# Solve
drop = ODEdrop2D(ic=(1,-1),t_end=600,het=g,V=V)
drop.solve()

# Make plots
fig = plt.figure(figsize=(6,3),tight_layout=True)
gs = GridSpec(2, 2, figure=fig)

ax_a = fig.add_subplot(gs[0],adjustable='box')
ax_b = fig.add_subplot(gs[1],adjustable='box')
ax_c = fig.add_subplot(gs[2],adjustable='box')
ax_d = fig.add_subplot(gs[3],adjustable='box')

# Part (a)
for tH in L:
    ax_a.plot(tH[0],tH[1],'k',lw=0.25)

plt.sca(ax_a)    
drop.draw(t_plots,ls='--',color='tab:blue')    
ax_a.set_xlim((-1.5,1.5))
ax_a.set_ylim([-0.05,1.2])
ax_a.set_ylabel('$h$')
ax_a.text(-2.4,ax_a.get_ylim()[1],'(a)')



# Part (b)
plt.sca(ax_b)
ax_b.plot(t_pde[::20],d_pde[::20],'k',lw=0.25,label='PDE')
ax_b.set_xlim((0,600))
ax_b.set_ylabel('$d$')
drop.plot(t=t_pde[::20],kind='width',ls='--',lw=1,label='ODEs')
ax_b.legend()
ax_b.text(-150,ax_b.get_ylim()[1],'(b)')

# Part (c)
x = np.linspace(-1.5,1.5,num=300)
ax_c.plot(x,g(x),'r',lw=0.5)
ax_c.set_xlim((-1.5,1.5))
ax_c.set_ylabel('$\\theta$')
ax_c.set_xlabel('$x$')
ax_c.text(-2.4,ax_c.get_ylim()[1],'(c)')


# Part (d)
plt.sca(ax_d)
ax_d.plot(t_pde[::20],theta_pde[::20],'k',lw=0.25)
ax_d.set_xlim((0,600))
ax_d.set_ylabel('$\\vartheta$')
ax_d.set_xlabel('$t$')
drop.plot(t=t_pde[::20],kind='angle',ls='--',lw=1)
ax_d.text(-150,ax_d.get_ylim()[1],'(d)')

plt.savefig('Figure03.png',bbox_inches='tight', dpi=200)