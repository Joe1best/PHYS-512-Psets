import numpy as np
import matplotlib.pyplot as plt
import random as rn
import matplotlib.animation as animation

import particle as P
import NBody as nb

def animate(i):
    global g,ax,fig,energy_file,pos_file,N
    g.evolve(nsteps=10,file_save=energy_file,file_save_pos=[pos_file,N_track])
    ptcl.set_data(g.densities)
    ptcl.set_cmap(plt.get_cmap('inferno'))
    return ptcl,

boundary_type ='Non-Periodic'
npart = 2**17
gridsize = 2**9
size = (gridsize,gridsize)
velInit = 0 
mass = 1/npart
init_mas = [mass for t in range(npart)]
soft = 0.8
s = P.system_init(npart,size,init_mas,npart_specific=None,npart_specificVel=velInit,boundary_type=boundary_type)

energy_file = open('Part3_NonPeriodic_Energy.txt','w')
pos_file = open('Part3_NonPeriodic_Position.txt','w')
N_track = 40


dt = 5
g = nb.NBody(size,s,dt,soft=soft,boundary_type=boundary_type)

niter = 300

Title = 'Part3_NonPeriodic.gif'
T = r'Simulation with $2^{17}$ Particles with $m$=7.6e-6, $dt$=5 with 300 frames'
x = []
y = []

labelsize = 15
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111, autoscale_on=False, xlim = (0,size[0]),ylim=(0,size[0]))
ax.tick_params(labelsize=labelsize)
ax.set_xlabel("X Position",fontsize=labelsize)
ax.set_ylabel("Y Position",fontsize = labelsize)
ax.set_title(T,fontsize=labelsize)

ptcl = ax.imshow(g.densities,origin='lower',vmin=g.densities.min(),vmax=g.densities.max())
plt.colorbar(ptcl)


an = animation.FuncAnimation(fig,animate,frames=niter,interval=10,repeat=True)
an.save(Title, writer='imagemagick')