import numpy as np 
import matplotlib.pyplot as plt
import time

def cylinder(n,R,pot):
    """
    Function that initializes the boundary conditions and the mask for a 3-D 
    cylinder of radius R and lenght L. This is done so the cylinder is centered
    in  a nxnxn space. 
    Inputs: 
        - n (int): resolution of the box
        - R (float): radius of the cylinder
        - pot (float): potential of the surface of the cylinder
    Outputs: 
        - V (array): boundary condition of the defined space along with the cylinder
        - bc (array): boundary condition
        - mask (array): mask of the cylinder
        - x, y (array): defines the grid of our space
    """
    midSlice = n//2
    s = np.linspace(0,n,n)
    x,y = np.meshgrid(s,s)

    V = np.zeros(x.shape)
    bc = 0*V 

    mask = np.zeros([n,n],dtype='bool')
    mask[0,:] = True
    mask[-1,:] = True
    mask[:,0] = True
    mask[:,-1] = True

    center_x, center_y = n//2,n//2
    cond = (x-center_x)**2+(y-center_y)**2 <= R**2
    mask[cond] = True
    bc[cond] = pot
    V = bc.copy()
    
    return V,bc,mask,x,y

def get_rhs(bc,mask,copy=False):
    """
    Makes the right hand side of the matrix formulation. 
    """
    rhs=np.zeros(bc.shape)
    
    rhs[:,:-1]= rhs[:,:-1] + bc[:,1:]
    rhs[:,1:]= rhs[:,1:] + bc[:,:-1]
    rhs[:-1,:]= rhs[:-1,:] + bc[1:,:]
    rhs[1:,:]= rhs[1:,:] + bc[:-1,:]
    rhs[mask]=0  
    
    return rhs 

def get_laplacian(A,mask,copy=False):
    """
    Calculates the laplacian operator.
    Inputs:
        - A (array): potential array 
        - mask (array): Array of bools where True is where the bc apply
    Ouputs: 
        - temp (array): laplacian operator in its matrix form 
    """
    if copy:
        A = A.copy()
    A[mask] = 0  # drop boundaries
    
    temp = 4*A  
    temp[:, :-1] -= A[:, 1:]  
    temp[:, 1:] -= A[:, :-1]  
    temp[:-1, :] -= A[1:, :]  
    temp[1:, :] -= A[:-1, :]  
    temp[mask] = 0  
    
    return temp

def true_V(x,y,R,n,pot,bc,mask):
    """
    Calculates the real potential using the expected formula of the potential
    """
    """
    cx,cy = n//2,n//2
    cond = (x-cx)**2+(y-cy)**2 <= R**2

    X,Y = x.ravel()-cx, y.ravel()-cy
    norm = np.sqrt(X**2+Y**2)-R
    true_V = np.log(norm)
    true_V[true_V<0] = 0
    true_V = true_V.reshape(x.shape)
    temp_V = true_V
    slope=1/temp_V[cx,0]
    true_V = -slope*true_V+pot
    true_V[cond] = pot
    """
    dv = bc[0, 0] - pot
    logr = np.log(np.sqrt(x**2 + y**2))
    lam = dv / (logr[0, 0] - logr[n//2, n//2+R])
    const = pot - lam*logr[n//2, n//2+R]
    true_V = lam*logr + const
    true_V[mask] = bc[mask]

    return true_V

def density(V,mask):
    """
    Computes the charge density associated with this potential
    Inputs: 
        - V (array): The solved potential of the grid
    Ouputs: 
        - density (array): The charge density on the grid
    """
    return get_laplacian(V,mask,copy=True)[1:-1,1:-1]


def relaxation(V,bc,mask,maxIter=10000,thresh=1e-2):
    """
    Implements the relaxation method as seen in class.
    Inputs:
        - V (array): potential grid that is solved
        - bc (array): boundary condition array
        - maxIter (int): maximum iterations of the relexation method
        - thresh (float): accuracy to which we want to solve in
    """
    st = time.time()
    r=get_rhs(bc,mask)-get_laplacian(V,mask,copy=True)
    rtr = np.sum(r*r)
    for i in range(maxIter):
        V[1:-1,1:-1] = (V[1:-1,:-2]+V[1:-1,2:]+V[2:,1:-1]+V[:-2,1:-1])/4.0
        V[mask] = bc[mask]
        r = get_rhs(bc,mask) - get_laplacian(V,mask,copy=True)
        rtr=np.sum(r*r)
        if rtr>thresh:
            continue
        elif rtr < thresh: 
            fi = time.time()-st
            print (f"The relaxation method converged after {i} iterations."
                  f"This algorithm converged in {fi}s")
            break
        elif i ==maxIter-1:
            fi = time.time()-st
            print (f"The relaxation did not converge before reaching the maximum"
                  f"iteration. Either increase the steps or change the intial bc's."
                  f"This took about {fi}s to finish. Potential may be solved poorly")
            break

def plot_three_results(V,V_true,d,sp,x,y,figsize=(20,6)):
    """
    Function that plots the numerical solutions, the analytic solutions and density in three
    side by side plots. 
    Inputs:
        - V (array): array  of solved potential 
        - V_true (array): theoretical expectations of the potential
        - d (array): density solved 
        - sp (int): sparness of data
        - x,y (array): the mesh of our space
    Outputs: 
        - graphs described above 
    """
    V_fiveP = V[::sp,::sp]
    E = np.gradient(V_fiveP)
    V_true_fiveP = V_true[::sp,::sp]
    E_true_fiveP = np.gradient(V_true_fiveP)
    
    fig, ax = plt.subplots(1,3,figsize=figsize)
    ax0 = ax[0].pcolormesh(V, vmin=0,vmax=1)
    ax[0].quiver(x[::10,::10].ravel(),y[::10,::10].ravel(),-E[1],-E[0])
    fig.colorbar(ax0, ax=ax[0])
    ax[0].set_title('Numerical Solution')
    ax1 = ax[1].pcolormesh(V_true,vmin=0,vmax=1)
    ax[1].quiver(x[::10,::10].ravel(),y[::10,::10].ravel(),-E_true_fiveP[1],-E_true_fiveP[0])
    ax[1].set_title('Analytic solution')
    fig.colorbar(ax1, ax=ax[1])
    ax2 = ax[2].pcolormesh(d,vmin=0,vmax=1)
    fig.colorbar(ax2, ax=ax[2])
    ax[2].set_title('Density')
    ax[0].set_xlabel("X pixels")
    ax[0].set_ylabel("Y pixels")
    ax[1].set_xlabel("X pixels")
    ax[2].set_xlabel("X pixels")
    plt.show()