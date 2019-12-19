import numpy as np
import matplotlib.pyplot as plt
import random as rn
import matplotlib.animation as animation

import particle as P
import NBody as nb

npart = 2**10
gridsize = 2**5
size = (gridsize,gridsize)
velInit = 0 
mass = 1
init_mas = [mass for t in range(npart)]
soft = 0.8
dt = 10

s = P.system_init(npart,size,init_mas,npart_specific=None,npart_specificVel=velInit)
g = nb.NBody(size,s,dt,soft=soft)
densities = g.densities
pos = g.posP
x = np.arange(0,gridsize)
y = np.arange(0,gridsize)
plt.imshow(densities,vmin=0,vmax=densities.max())
plt.plot(pos[:,1],pos[:,0],'.',color='red')
for i in range(gridsize):
    plt.axvline(x[i],color='black')
    plt.axhline(y[i],color='black')
plt.colorbar()
plt.show()