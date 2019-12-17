import numpy as np
import matplotlib.pyplot as plt
import random as rn
import matplotlib.animation as animation

import particle as P
import NBody as nb

def animate(i):
    global system,ax,fig
    g.evolve()
    ptcl.set_data(g.posP[:,0],g.posP[:,1])
    return ptcl,

npart = 1
gridsize = 2**9
size = (gridsize,gridsize)

mass = 10
init_mas = [mass for t in range(npart)]

npart_specific = np.array([[gridsize/2,gridsize/2]])
npart_specificVel = np.array([[0,0]])

s = P.system_init(npart,size,init_mas,npart_specific=npart_specific,npart_specificVel=npart_specificVel)

dt = 1
g = nb.NBody(size,s,dt,soft=0.1)

niter = 100

Title = 'Part1.gif'
T = f'Stationary object with $v_o$=0, $m$={mass}, $dt$={dt} with {niter} frames'
x = []
y = []

labelsize = 15
fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111, autoscale_on=False, xlim = (0,size[0]),ylim=(0,size[0]))
ax.tick_params(labelsize=labelsize)
ax.set_xlabel("X Position",fontsize=labelsize)
ax.set_ylabel("Y Position",fontsize = labelsize)
ax.set_title(T,fontsize=labelsize)

ptcl, = ax.plot([],[],'*',markersize=10,color='black')

an = animation.FuncAnimation(fig,animate,frames=niter,interval=dt,repeat=True)

an.save(Title, writer='imagemagick')