# -*- coding: utf-8 -*-
"""
Data generation for Figure 05
"""
import sys
sys.path.append("../main/")
from ODEdrop2D import *
from pdeloader import *
from matplotlib.gridspec import GridSpec
from scipy.io import loadmat
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle

V = lambda t: Vperiodic(t,m=20,p=300,Vavg=4,Vamp=0.75)
g = lambda x: 1 + .2*np.tanh(20*np.cos(8*np.pi*x))

# Solve the ODE problem
drop = ODEdrop2D(ic=(1.,-1),t_end=600,het=g,V=V)
drop.solve()

# Load PDE data
t_pde, d_pde, l_pde = pdeloader("Figure05_PDE.mat")

fig = plt.figure(figsize=(8,2.5))
gs = GridSpec(1, 2, figure=fig,wspace=0.24)

# Plot profiles
ax_a = fig.add_subplot(gs[0],adjustable='box')
plt.sca(ax_a)
drop.draw(np.arange(225,400,25),xlim=(0,2.5),color='tab:blue')

for ui,ri,si in zip([0.8,1.2],[0.3,0.6],['$\\vartheta_{min}$','$\\vartheta_{max}$']):
    u = np.linspace(np.pi,np.pi-np.arctan(ui))
    ax_a.plot(2.44+ri*np.cos(u),ri*np.sin(u),'k')
    ax_a.plot([2.44+(ri+0.1)*np.cos(u[-1]),2.44],[(ri+0.1)*np.sin(u[-1]),0],'k')
    ax_a.text(2.44+ri*np.cos(u[20]),ri*np.sin(u[20]),si,horizontalalignment='right')
ax_a.set_xlim((0,2.5))
ax_a.set_xlabel('$x$')
ax_a.set_ylabel('$h$')
ax_a.text(-0.5,ax_a.get_ylim()[1],'(a)')

# Plot evolution of the radius
ax_b = fig.add_subplot(gs[1],adjustable='box')
ax_b.set(xlim=(0,500),ylim=(1,2.5))

t = np.arange(501)
sol = drop.evaluate(t)
ax_b.plot(t_pde[::20],d_pde[::20],'k',lw=0.25,label='PDE')
ax_b.plot(t,0.5*(sol[0]-sol[1]),'--',lw=1,label='ODEs')
ax_b.legend()
ax_b.add_patch(Rectangle((150, 2.430),300, .015, ls="-", lw=0.25,ec="r"))
ax_b.set_xlabel('$t$')
ax_b.set_ylabel('$d$')
ax_b.text(-110,ax_b.get_ylim()[1],'(b)')

plt.plot([163,150,np.nan,463,450],[2.25,2.430,np.nan,2.25,2.430],':')

axins = inset_axes(ax_b, width="100%", height="100%",
                   bbox_to_anchor=(.3, .3, .6, .5),
                   bbox_transform=ax_b.transAxes,loc=3)
axins.plot(t,0.5*(sol[0]-sol[1]),'--',lw=1,label='ODEs')
axins.plot(t_pde[::20],d_pde[::20],'k',lw=0.25,label='PDE')
axins.set(xlim=(150,450),ylim=(2.430,2.445))


plt.savefig('Figure05.png', bbox_inches='tight',dpi=200)