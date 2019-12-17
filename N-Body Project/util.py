import numpy as np


def evolve(posP,velocityP,mass,f,dt,size):
    """
    Evolves the particle's position and momentum using the leap frog method seen in class
    Inputs: 
        - f_current (array): The forces on the particle 
        - f_new (array): changed forces on the particle 
        - dt (float): time step in seconds 
    """
    #print (velocity.shape,f.shape,mass.shape,"icitte")
    velocityP = velocityP+f*dt/mass

    # update position
    posP = posP+velocityP*dt
    posP = posP%size
    return posP,velocityP 

def loadPosition(file):
    """
    Function that loads in the position and converts them to an actual numpy 
    array 
    Input(s):
        - file (file): file containing the positions of the particle
    Output(s):
        - A (array): positions of the particles in a numpy array format
    """
    A = []
    for i in range(300):
        H = []
        for x in range(40):
            string = file.readline().replace('[[','').replace(']]','').replace('\n','').replace('[','').replace(']','')
            num = np.array(string.split())
            num = num.astype(np.float)
            H.append(num)
        A.append(H)
    return A

def calculate_dist(grid,pos):
    """
    Calculates the position between a grid point and a particle's position
    """
    return np.sqrt((grid[0]-pos[0])**2+(grid[1]-pos[1])**2)