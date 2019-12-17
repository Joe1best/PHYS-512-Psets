import util as ut
import numpy as np


class NBody: 
    def __init__(self,size,particleList,dt,soft=0.1,G=1,boundary_type='Periodic'):
        """
        The NBody class that specifies the simulation. 
        Input(s):
            - size (x,y): size of the grid 
            - particleList (system_init object): Initial list of the particles along with velocity and position
            - dt (float): step in time taken in the simulation 
            - soft (float): softner used for the green function (default 0.1)
            - G (float): Gravitational constant (default 1.0)
            - boundary type: Either set to Periodic or Non-Periodic
        """

        self.boundary_type = boundary_type

        #If the boundary type is non-periodic, we multiply the grid by 2 in order to make
        #the particles not feel any from the other side. Everything will only be generated 
        #in the orignal quadrant (all the particles live in the original (size[0],size[1]))
        if self.boundary_type == 'Non-Periodic':
            self.size = (2*size[0],2*size[1])   
        else:
            self.size = size
        self.soft = soft
        self.G = G
        self.dt = dt
        self.posP = particleList.pos
        self.velocityP = particleList.velocities
        self.mass = particleList.masses
        x, y = np.arange(self.size[0],dtype=float), np.arange(self.size[0],dtype=float)
        self.mesh = np.array(np.meshgrid(x,y))
        self.density_assignment()

        self.green()
        
       
    def density_assignment(self):
        """
        Function that assigns the density of the grid according to the 
        Nearest Grid Points (NGP) scheme. A trial with the CIC model (cloud-in-cell)
        scheme lead to a lot of problems with the inverse scheme, so I opted for the NGP
        This scheme is just done by assigning the mass of any particle to its nearest 
        gridpoint (since the cells are of unit length, no need to divide by anything)
        """
       #Assigns to nearest grid point
        self.intPos = np.rint(self.posP).astype('int') % self.size[0]
        
        #To make it easier to access the gridpoints later on 
        self.mesh_modified = tuple(self.intPos[:, i] for i in range(2))

        #Instead of looping (to save time), I found this method online to bin the masses 
        #in their nearest gridpoint
        E = np.linspace(0, self.size[0]-1, num=self.size[0]+1)
        E = np.repeat([E], 2, axis=0)
        hist = np.histogramdd(self.intPos, bins=E,weights=self.mass.flatten())
        
        self.densities = hist[0]
        
    def green(self):
        """
        Function that defines the greenfunction. 
        Note: this will only work for square grids as there was not enough time 
        to implement this with rectangular grids. 
        """
        r = np.sum(self.mesh**2,axis=0)
        r[r<self.soft**2] = self.soft**2
        r += self.soft**2
        r = np.sqrt(r)
        
        g = 1/(4*np.pi*r)
        
        #To get periodicity, we flip the corners to get the same behavior around
        h_x,h_y = self.size[0]//2, self.size[1]//2

        try:
            g[h_x:, :h_y] = np.flip(g[:h_x,:h_y],axis=0)
            g[:,h_y:] = np.flip(g[:,:h_y],axis=1)
        except: 
            g[h_x:, :h_y+1] = np.flip(g[:h_x+1,:h_y+1],axis=0)
            g[:,h_y:] = np.flip(g[:,:h_y+1],axis=1)
        self.g = g
                    
    def pot(self):
        """
        Function that retrieves the potential of the grid by convoluting the density 
        of the grid and the green function of the grid. 
        Made sure to take of the descrepancy between the fft and the center of the particles
        Without the for loops (with the np.roll), the particles were all attracted to themselves
        and started to drift towards the origin.
        """
        
        ffD = np.fft.rfftn(self.densities)
        ffG = np.fft.rfftn(self.g)
        ffV = ffD*ffG
        
        V = np.fft.irfftn(ffV)

        #Need to shift and average the potential to center it back to particle
        for i in range(2):
            V = 0.5*(np.roll(V,1,axis=i)+V)
        
        #If the boundary conditions is not periodic, set the potential to 0
        if self.boundary_type == 'Non-Periodic':
            V[0:,0] = 0
            V[0,-1] = 0
            V[-1:,-1] = 0
            V[-1:,0] = 0 
        return V 
    
    def forces_mesh(self): 
        """
        To get the forces, we simply take the gradient of the potential. 
        To take the gradient, we use the central difference
        """
        fmesh = np.zeros([2,self.size[0],self.size[1]])
        V = self.pot()
        fmesh[0] = 0.5*(np.roll(V,1,axis=0)-np.roll(V,-1,axis=0))
        fmesh[1] = 0.5*(np.roll(V,1,axis=1)-np.roll(V,-1,axis=1))
        
        #Multiply by gravitational constant and the densities to get the force
        fmesh = -fmesh*self.densities*self.G
        return fmesh 
    
    def forces_pctls(self):
        """
        Function to interpolate the forces using the inverse scheme of the 
        NGP density scheme. 
        """
        
        fxy = np.moveaxis(self.forces_mesh(),0,-1)

        f_final = fxy[self.mesh_modified]
         
        return f_final
    
    def totalEnergy(self):
        """
        Function to compute the total energy of the system using
        E = V+0.5*m*v**2
        """
        K = np.sum(self.mass*self.velocityP**2)
        P = -0.5*np.sum(np.sum(self.pot())*self.densities)
        T = K + P
        return T 
    
    def evolve(self,nsteps=1,file_save=None,file_save_pos=None):
        """
        Evolves the system and saves the energy along with position for further 
        analysis. The evolution uses the leapfrog method, which is written in the 
        util.py file
        Input(s):
            - nsteps (int): how many steps per evolution before it saves a result
            - file_save (file): File in which we want to save the total energy of the system
            - file_save_pos (file,int): First entry is the file in which we want to save the 
            energy is, second entry is the number of particles we want to track.
        """
        for i in range(nsteps):
            F = self.forces_pctls()
            self.posP,self.velocityP = ut.evolve(self.posP,self.velocityP,self.mass,F,self.dt,self.size[0])

            self.density_assignment()

        energy = self.totalEnergy()
        if file_save is not None:
            file_save.write(f"{energy}\n")
            file_save.flush()
        if file_save_pos is not None:
            file_save_pos[0].write(f"{self.posP[:file_save_pos[1]]}\n")
            file_save_pos[0].flush()
        return energy,self.posP