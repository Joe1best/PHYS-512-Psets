import numpy as np
import matplotlib.pyplot as plt
import random as rn
import matplotlib.animation as animation
from matplotlib.colors import LogNorm, Normalize


import particle as P
import NBody as nb

def animate(i):
    global g,ax,fig,energy_file
    g.evolve(nsteps=10,file_save=energy_file)
    ptcl.set_data(g.densities)
    ptcl.set_cmap(plt.get_cmap('inferno'))
    return ptcl,

npart = 2**17
gridsize = 2**9
size = (gridsize,gridsize)
velInit = 0 
mass = 40
init_mas = [mass for t in range(npart)]
soft = 10
s = P.system_init(npart,size,init_mas,npart_specific=None,npart_specificVel=velInit,soft=soft,cosmos=True)
energy_file = open('Part4.txt','w')
dt = 330
g = nb.NBody(size,s,dt,soft=soft)

niter = 450

Title = 'Part4.gif'
T = r'Simulation with $2^{17}$ Particles with Mass Fluctuations $\propto$ to $k^{-3}$, $dt$=330 with 450 frames'
x = []
y = []

labelsize = 15
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111, autoscale_on=False, xlim = (0,size[0]),ylim=(0,size[0]))
ax.tick_params(labelsize=labelsize)
ax.set_xlabel("X Position",fontsize=labelsize)
ax.set_ylabel("Y Position",fontsize = labelsize)
ax.set_title(T,fontsize=labelsize)

gCopy = g.densities.copy()
gCopy[gCopy==0] = gCopy[gCopy!=0].min()*1e-3

ptcl = ax.imshow(g.densities,origin='lower',norm=LogNorm(vmin=gCopy.min(),vmax=gCopy.max()))
plt.colorbar(ptcl)

an = animation.FuncAnimation(fig,animate,frames=niter,interval=10,repeat=True)

an.save(Title, writer='imagemagick')