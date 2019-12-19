# N-Body Simulation 

This is the repository for the N-Body simulation done for PHYS-512 final project. In this project, we were tasked to answer questions 
according to the [Project Guidelines](https://github.com/Joe1best/PHYS-512-Psets/blob/master/N-Body%20Project/project_guidelines.pdf). 

To approach this project, I used the Particle Mesh (PM) method described in more details below. Important to note that this project is only
implemented in 2-D for the time being. Further work will see this implemented in N-dimensions. The results for each question can be found in the [gifs folder](https://github.com/Joe1best/PHYS-512-Psets/tree/master/N-Body%20Project/gifs). Methods of setting up and running the simulation can be found in the python file allocated for each part.

## Particle Mesh 

Below, I describe the steps taken to accomplish the simulation. The crux of this method is to numerically integrate Poisson's equation; 

<a align ="center" href="https://www.codecogs.com/eqnedit.php?latex=\Delta&space;\Phi(\vec{g}_{k,l,m})=4\pi&space;G\rho(\vev{g}_{k,l,m})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Delta&space;\Phi(\vec{g}_{k,l,m})=4\pi&space;G\rho(\vev{g}_{k,l,m})" title="\Delta \Phi(\vec{g}_{k,l,m})=4\pi G\rho(\vev{g}_{k,l,m})" /></a>

<a align = "center" href="https://www.codecogs.com/eqnedit.php?latex=\vec{F}(\vec{g}_{k,l,m})&space;=&space;-m&space;\nabla\Phi(\vec{g}_{k,l,m})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\vec{F}(\vec{g}_{k,l,m})&space;=&space;-m&space;\nabla\Phi(\vec{g}_{k,l,m})" title="\vec{F}(\vec{g}_{k,l,m}) = -m \nabla\Phi(\vec{g}_{k,l,m})" /></a>

To do this, we must;

1) Calculate the mass density on grid 

2) Solve Poisson's equation on grid 

3) Differentiate Potential to get forces 

4) Interpolate forces back to particles

This sounds like a waste of time and computer resources. However, this is exceptionally fast in practice. To implement this, all the functions necessary are located in [NBody.py](https://github.com/Joe1best/PHYS-512-Psets/blob/master/N-Body%20Project/NBody.py)

### Calculating mass density on grid 

There are multiple schemes out there to calculate the grid density. Some of which include the Nearest Grid Point (NGP), Cloud-In-Cell
(CIC), or Triangular-Shaped-Cloud (TSC). The chosen sheme in this project was the NGP one. Although this method is too crude and generally not 
smooth, it was easy and fast to implement and requires less running time. The CIC method is a common choice for N-Body simulation, however, it takes more time to run for 
these types of simulations. The TSC scheme is the smoothest but hardest to implement. I will go over the NGP scheme down below.

Consider a particle that has a position (x,y,z) in a grid that has spacing H. This particle's mass will be binned in the closest grid point defined by <a href="https://www.codecogs.com/eqnedit.php?latex=\vec{g}_{k,l,m}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\vec{g}_{k,l,m}" title="\vec{g}_{k,l,m}" /></a>. Therefore, each gridpoint will have a total mass that is equal to the sum of their binned particles. We then divide by H to get the density of each gridpoint. In this project, H was unitary (defined to be equal to 1). Hence the density as a function of gridpoints was just the sum of the binned masses. The function that takes care of this is called **density_assignment()** found [here](https://github.com/Joe1best/PHYS-512-Psets/blob/master/N-Body%20Project/NBody.py#L40). Below is the result of this function with 1024 particles each with mass 1 in a 32x32 grid, 

![Density](https://github.com/Joe1best/PHYS-512-Psets/blob/master/N-Body%20Project/Documentation/Figure_1.png)

As one can see, the density is not very smooth when using the NGP scheme. However, as we increase the size of the grids, the behavior becomes smooth enough to approximately simulate correctly the interactions between those particles. 

### Solve Poisson's equation on grid 

We want to solve the following equation for the potential,

<a href="https://www.codecogs.com/eqnedit.php?latex=\Delta&space;\Phi_{k,l,m}&space;=&space;\rho_{k,l,m}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Delta&space;\Phi_{k,l,m}&space;=&space;\rho_{k,l,m}" title="\Delta \Phi_{k,l,m} = \rho_{k,l,m}" /></a>

This will be done using the Green's function method. Green's function on the grid can be proven to be equal to, 

<a href="https://www.codecogs.com/eqnedit.php?latex=G(\vec{g}_{k,l,m})&space;=&space;\frac{1}{4\pi\vec{r}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?G(\vec{g}_{k,l,m})&space;=&space;\frac{1}{4\pi\vec{r}}" title="G(\vec{g}_{k,l,m}) = \frac{1}{4\pi\vec{r}}" /></a> 

where <a href="https://www.codecogs.com/eqnedit.php?latex=\vec{r}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\vec{r}" title="\vec{r}" /></a> is the distance from the origin to the grid point. Important to note however, that we will be exploring periodic boundary conditions and non-periodic boundary conditions. For Non-Periodic, green's function looks as one can expect (left in the figure below), 

![Green](https://github.com/Joe1best/PHYS-512-Psets/blob/master/N-Body%20Project/Documentation/Green_both.PNG)

However, for Periodic behavior, we need to flip the function around each corner to make it look like the right figure above. Finally, we add a softening constant that makes the green function continuous and not having any spike behavior (usually around the zero mark). For the above figure, the softening constant used was 10, although the norm used in this project is about 0.8 to 1.0 . All of this is taken care of in the **green()** function found [here](https://github.com/Joe1best/PHYS-512-Psets/blob/master/N-Body%20Project/NBody.py#L62).

Having the green function, the potential is just equal to the convolution of the density function found previously with the green function. In Fourier Space, this just represents a simple multiplication between the Fast Fourier tranfrom of the density of the grid multiplied by the Fast Fourier Transfrom of the green function on the grid, the result of which will be inversed Fourier Transfromed. This is done in the **pot()** function [here](https://github.com/Joe1best/PHYS-512-Psets/blob/master/N-Body%20Project/NBody.py#L86).

### Differentiating potential to get forces

Since the force is equal to the negative gradient of the function, we take the gradient of the found potential using what is the called the central difference method defined below, 

<a href="https://www.codecogs.com/eqnedit.php?latex=F_x(\vec{g}_{k,l,m})=-\frac{1}{2}G_o\rho(\vec{g}_{k,l,m})(\Phi(\vec{g}_{k&plus;1,l,m})-\Phi(\vec{g}_{k-1,l,m}))\\&space;F_y(\vec{g}_{k,l,m})=-\frac{1}{2}G_o\rho(\vec{g}_{k,l,m})(\Phi(\vec{g}_{k,l&plus;1,m})-\Phi(\vec{g}_{k,l-1,m}))\\&space;F_z(\vec{g}_{k,l,m})=-\frac{1}{2}G_o\rho(\vec{g}_{k,l,m})(\Phi(\vec{g}_{k,l,m&plus;1})-\Phi(\vec{g}_{k,l,m-1}))\\" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F_x(\vec{g}_{k,l,m})=-\frac{1}{2}G_o\rho(\vec{g}_{k,l,m})(\Phi(\vec{g}_{k&plus;1,l,m})-\Phi(\vec{g}_{k-1,l,m}))\\&space;F_y(\vec{g}_{k,l,m})=-\frac{1}{2}G_o\rho(\vec{g}_{k,l,m})(\Phi(\vec{g}_{k,l&plus;1,m})-\Phi(\vec{g}_{k,l-1,m}))\\&space;F_z(\vec{g}_{k,l,m})=-\frac{1}{2}G_o\rho(\vec{g}_{k,l,m})(\Phi(\vec{g}_{k,l,m&plus;1})-\Phi(\vec{g}_{k,l,m-1}))\\" title="F_x(\vec{g}_{k,l,m})=-\frac{1}{2}G_o\rho(\vec{g}_{k,l,m})(\Phi(\vec{g}_{k+1,l,m})-\Phi(\vec{g}_{k-1,l,m}))\\ F_y(\vec{g}_{k,l,m})=-\frac{1}{2}G_o\rho(\vec{g}_{k,l,m})(\Phi(\vec{g}_{k,l+1,m})-\Phi(\vec{g}_{k,l-1,m}))\\ F_z(\vec{g}_{k,l,m})=-\frac{1}{2}G_o\rho(\vec{g}_{k,l,m})(\Phi(\vec{g}_{k,l,m+1})-\Phi(\vec{g}_{k,l,m-1}))\\" /></a>

where <a href="https://www.codecogs.com/eqnedit.php?latex=G_o" target="_blank"><img src="https://latex.codecogs.com/gif.latex?G_o" title="G_o" /></a> is the Gravitational constant (will be defined as 1 in this project). Since only 2-D was done, the first two equations were used without the m. This is implemented in **forces_mesh()** [here](https://github.com/Joe1best/PHYS-512-Psets/blob/master/N-Body%20Project/NBody.py#L113). 

### Interpolate forces back to particles 

Easiest step of them all. It just involves doing the inverse scheme of the density assignment. In the case of the NGP, it just involves extending the forces felt by an arbitrary gridpoint to all the particles binned in that gridpoint. The function **forces_ptcl()** [here](https://github.com/Joe1best/PHYS-512-Psets/blob/master/N-Body%20Project/NBody.py#L127) takes care of that. 


