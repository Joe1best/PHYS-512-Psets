import random as rn
import numpy as np
import matplotlib.pyplot as plt



class particle: 
    def __init__(self,m,x,y,vx=0,vy=0): 
        """
        Defines the particle. Each particle has a mass, a position and momentum
        Inputs: 
            - m (float): mass of the particle 
            - x,y (float): x and y components respectively of the particle's position 
            - vx,vy (float): x and y components respectively of the particle's velocity
        """
        self.mass = m 
        self.position = (x,y)
        self.velocity = (vx,vy)

class system_init: 
    def __init__(self,npart,size,init_mass,npart_specific=None,npart_specificVel=None,boundary_type='Periodic',soft=None,cosmos=False):
        """
        Defines the total system of particles and puts them in a list. 
        Inputs: 
            - npast (int): total number of particles 
            - size (int): (x,y) of the size of the grid 
            - init_mass (array): initial mass of the particles
            - npart_specific (array): inital position of the particles (specified)
            - npart_specificVel (array): inital velocity of the particles (specified)
            If the last two are not specified, then they are generated below via the two 
            functions
        """
        self.boundary_type = boundary_type
        def initial_pos(size,npart,npart_specific=None): 
            """
            Depending on the grid size, generates random position on the grid. 
            Inputs:
                - size (array): size of grid. In the format of (x,y)
                - npart (int): total number of particles of the system
                - npar_specific (array): specific boundary condition pre-defined by the user
            """
            self.cosmos = cosmos
            if self.boundary_type == 'Periodic':
                init_cond = []
                if npart_specific is None: 
                    l = 0
                else: 
                    l = len(npart_specific)
                    for k in range(l):
                        pos = (npart_specific[k][0], npart_specific[k][1])
                        init_cond.append(pos)

                for p in range(npart-l):
                    pos = (rn.random())*(size[0]-1), (rn.random())*(size[1]-1)
                    init_cond.append(pos)
            elif self.boundary_type == 'Non-Periodic':
                xmin,xmax = 1,size[0]-1
                init_cond = []
                if npart_specific is None: 
                    l = 0
                else: 
                    l = len(npart_specific)
                    for k in range(l):
                        if npart_specific[k][0].max()>xmax or npart_specific[k][1].max()>xmax or npart_specific[k][0].min()<xmin or npart_specific[k][1]<xmin: 
                            raise ValueError(f'The position of the Particle must be within the boundary of 1 to {size[0]-1}')
                        
                        pos = (npart_specific[k][0], npart_specific[k][1])
                        init_cond.append(pos)
                for p in range(npart-l):
                    pos = rn.uniform(1.0001,size[0]-1.0001), rn.uniform(1.0001,size[1]-1.0001)
                    init_cond.append(pos)
            return init_cond
        
        def initial_vel(size,npart,maxSpeed,npart_specific=None):
            """
            If the velocities are not specified (using npart_specific), this function will
            generate random velocity according to a Gaussian with a std of 1. 
            If npart_specific is given to be 0, all the particles will have an inital velocity of 0 
            Input(s):
                - size (x,y): size of the grid 
                - npart (int): the number of particles 
                - maxSpeed (float): Maximum intial speed of the particles
                - npart_specific (array): Specified inital conditions of the particles
            """
            init_cond = []
            
            if np.isscalar(npart_specific)==False:
                if npart_specific is None: 
                    l = 0
                else: 
                    l = len(npart_specific[0])
                    for k in range(l):
                        vel = (npart_specific[k][0], npart_specific[k][1])
                        init_cond.append(vel)
            
                for p in range(npart-l):
                    vel = (np.random.normal()*maxSpeed, np.random.normal()*maxSpeed)
                    init_cond.append(vel)
            else: 
                for p in range(npart):
                    vel = (0,0)
                    init_cond.append(vel)
            return init_cond
        
        def initial_mass(size,init_mass,posP=None,soft=None):
            if self.cosmos == False:
                return np.array([init_mass.copy()]).T
            else: 
                posP = np.rint(posP).astype('int') % size[0]
                
                #To find the power spectrum
                k_x = np.real(np.fft.fft(posP[:,0]))
                k_y = np.real(np.fft.fft(posP[:,1]))
                k_total = np.sqrt(k_x**2+k_y**2)
                k_total[k_total<soft] = soft
                m = init_mass/k_total**3
                
                return np.array([m.copy()]).T


        self.nparticles = npart
        init_cond = initial_pos(size,npart,npart_specific=npart_specific)
        init_vel = initial_vel(size,npart,1,npart_specific=npart_specificVel)
        self.particles = np.asarray([particle(m,x[0],x[1],vx=v[0],vy=v[1]) for m,x,v in zip(init_mass,init_cond,init_vel)])   
        self.velocities = np.asarray([self.particles[i].velocity for i in range(npart)])
        self.pos = np.asarray([self.particles[i].position for i in range(npart)])
        self.masses = initial_mass(size,init_mass,posP=self.pos,soft=soft)