import numpy as np 
import matplotlib.pyplot as plt


def heat_solver(dx,dt,t_max,x_max,k,C):
    m = k*dt/dx**2
    if m>0.5:
        print ("m must be greater than 1/2 for the pde to converge.")
        exit(1)
    x, t = np.arange(0,x_max+dx,dx), np.arange(0,t_max+dt,dt)
    r, c = len(t), len(x)
    
    T = np.zeros([r,c])
    for i in range(0,r-1):
        for j in range(0,c-1):
            T[i,0] = t[i]*C
            T[i+1,j] = T[i,j] + m*(T[i,j-1] - 2*T[i,j] + T[i,j+1]) 

    return x,T,t    

def plot(x,T,t,dt,sp):
    colors = plt.cm.jet(np.linspace(0,1,len(t[::sp])))
    i = 0
    fig, ax = plt.subplots(figsize=(20,10))
    for times in t[::sp]: 
        index = int(times/dt)
        ax.plot(x,T[index,:],label= f't={times}',color=colors[i],linewidth=4)
        i+=1
    ax.set_xlabel("Position (m)",fontsize=15)
    ax.set_ylabel("Temperature (K)",fontsize=15)
    ax.set_title("Evolution of temperature along center of the box as a function of time",fontsize=15)
    ax.tick_params(labelsize=15)
    plt.legend(fontsize=15)
    plt.show()