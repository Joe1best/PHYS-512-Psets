import time 
import numpy as np 
import Q1

def cg(bc,mask,maxIter=10000,V0=None,thresh=1e-2):
    """
    this runs a conjugate gradient solver to solve Ax=b where A
    is the Laplacian operator interpreted as a matrix, and b is the contribution from the 
    boundary conditions.  Incidentally, one could add charge into the region by adding it
    to b (the right-hand side or rhs variable)
    """
    if V0 is None:
        V0 = bc.copy()
    st = time.time()
    r = Q1.get_rhs(bc,mask)-Q1.get_laplacian(V0,mask)
    p = r.copy()
    V = V0.copy()
    rtr = np.sum(r*r)
    for k in range(maxIter):
        Ap=Q1.get_laplacian(p,mask)
        alpha=  rtr/np.sum(Ap*p)
        V += alpha*p
        r -= alpha*Ap
        rtr_new = np.sum(r*r)
        beta = rtr_new/rtr
        p = r + beta*p
        rtr = rtr_new
        if rtr < thresh:
            fi = time.time()-st
            print (f"The conjugate gradient converged in {fi}s after {k} iterations")
            break
        elif k == maxIter-1:
            fi = time.time()-st
            print (f"The maximum number of iterations was reached before CG can converge."
                  f" Please increase the max number of iterations. The potential might be solved"
                  f" poorly. This took {fi}s")
            break 
    V[mask] = bc[mask] 
    return V, k