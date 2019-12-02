import time
from scipy.interpolate import griddata
import numpy as np
import Q2

"""
I would like to give credit where its due and thank Thomas Vandal for helping 
me get through this problem. It took me so long by myself and it was not working,
so he came to the rescue and gave me a huge hand with it. 
"""

def lower_res(space):
    """
    Lower the resolution by a factor of 2 along each direction. 
    NOTE: everytime this is called, space must have an even dimension 
    and should be equal in size. This means that we have to choose an 
    initial n such that it is always a multiple of 2, aka 2^(k).
    Inputs: 
        - space (mask or bc): initial mask or bc array 
    Outputs: 
        - new_space (mask or bc): initial mask or bc array but half the resolution
    """
    s = space.shape

    n = s[0]
    new_space = space.reshape(n//2, 2,-1,2).swapaxes(1,2).reshape(-1,2,2)
    m = np.max(np.abs(new_space),axis=(1,2))
    m_t = np.max(new_space,axis=(1,2))
    min_t = np.min(new_space,axis=(1,2))
    new_space = np.where(np.equal(m,m_t),m_t,min_t)
    new_space = new_space.reshape(n//2,n//2)
    
    return new_space

def higher_res(space):
    """
    Increase the resolution by a factor of 2 along each direction. Uses 
    griddata from scipy which interpolates unstructured D-dimensional data.
    Inputs:
        - space (mask or bc): initial mask or bc array 
    Ouputs: 
        - new_space (mask or bc): initial mask or bc array but double the resolution
    
    """
    s = space.shape
    arr = np.linspace(0,s[0]-1,num=s[0])
    x,y = np.meshgrid(arr,arr)
    p = np.asarray([x.ravel(),y.ravel()]).T
    
    arr = np.linspace(0,s[0]-1,num=2*s[0])
    x,y = np.meshgrid(arr,arr)
    p_2 = np.asarray([x.ravel(),y.ravel()]).T
    new_space = griddata(p,space.ravel(),p_2,method='cubic')
    new_space = new_space.reshape(2*s[0],2*s[0])
    
    return new_space

def res_cg(bc,mask,V0=None,npass=6,thresh=1e-2,maxIter=10000):
    """
    Implements conjugate gradient with changing resolution as seen in class.
    """
    st = time.time()
    all_mask = [None]*npass
    all_bc = [None]*npass
    all_V = [None]*npass
    all_mask[0] = mask
    all_bc[0] = bc
    
    for i in range(1,npass):
        all_mask[i]= lower_res(all_mask[i-1])
        all_bc[i]= lower_res(all_bc[i-1])
    
    all_V[-1], count = Q2.cg(all_bc[-1], all_mask[-1], V0=V0, maxIter=maxIter, thresh=thresh)

    for i in range(npass-2, -1, -1): #dont know why we do this, it was defined in class
        V0 = higher_res(all_V[i+1])
        all_V[i], k = Q2.cg(all_bc[i], all_mask[i], V0=V0, maxIter=maxIter, thresh=thresh)
        count += k
    fi = time.time() - st
    
    print(f"Changing the resolution, the conjugate gradient took {fi}s with {count} iterations")    
    V = all_V[0]  
    return V