# -*- coding: utf-8 -*-
"""
Created on Wednesday Sep 18 16:50:55 2019

@author: joebe
"""
import numpy as np


#This block is a script to run simpson's rule, but on half the interval
def simpson(f,a,b,fa,fb):
    """
    Function that integrates f the interval a to b by splitting it in n subinterval and using Simpson's rule 
        -Inputs: f (function) this is the function to be integrated 
                 a (float) Lower bound of the intergral
                 b (float) Upper bound of the integral 
                 fa (float) Value of f at a
                 fb (float) Value of f at b
        -Ouputs: I (float) Value of the integral between the inteval a to b
                 h (float) Midpoint of the function
                 fh (float) Value of f at h
    """
    #Juts checks a bunch of cases for the interval a to b
    if a == b:
        return 0
    if a < b:
        sign = 1
    else:
        sign = -1
    h = (a+b)/2    #Calculates the height 
    fh = f(h)      #The value at that height
    dx = (b-a)/(6) #usually the denominator is 3*n, but here n=2, since we are halving the interval
    I = dx*sign*(fa+4*fh+fb) #--> Simpson's rule 
    return I, h, fh 

def dynamic_simpson(f,a,b,h,fa,fb,fh,guess,tolerance,nmax,count=0):
    """
    Function that estimates the error of the integral using |S(a,h)+S(h,b)-S(a,b)|<15*tolerance (this is 
    from the wikipedia page about adaptative simpson's rule), where [a,b] is our integral, and S is 
    the value of simspon rule. 
    If this error is smaller than the tolerance, it gives back the normal value of the intgral using 
    simpson(). If not, it splits the interval into another half and tests the error again.
        -Inputs: f (function) Function that is integrated 
                 a (float) Lower bound of the integral 
                 b (float) Upper bound of the integral 
                 h (float) Half point of the interval 
                 fa (float) Value of f at a
                 fb (float) Value of f at b
                 fh (float) Value of f at h 
                 guess (float) Given by Simpson's rule on the whole interval 
                 tolerance (float) Tolerance on how precise the integral is
                 nmax (int) maximum number of interation before exiting the recursion
                 count (int) number of iterations that it took for the integral to converge 
    
        -Ouputs: I (float) integral value using simpson's rule
                 Iter (int) how many iterations it took for the intergal to converge using the recursive
                 algorithm below. 
    """
    if count < nmax:
        Guess_ah,l_h,l_fh = simpson(f,a,h,fa,fh)
        Guess_hb,r_h, r_fh = simpson(f,h,b,fh,fb)
        error = np.abs(Guess_ah+Guess_hb - guess)
        if error < 15*tolerance:
            return guess, count
        else:
            guess_l, count_l = dynamic_simpson(f,a,h,l_h,fa,fh,l_fh,Guess_ah,tolerance/2,nmax,count+1)
            guess_r, count_r = dynamic_simpson(f,h,b,r_h,fh,fb,r_fh,Guess_hb,tolerance/2,nmax, count+1)
            
            guess_t = guess_l+guess_r
            count_t = count_l+count_r
            return guess_t,count_t
    else:
        print(f'Sadly my integral did not converge in {nmax} steps...oh well :(. Maybe increase nmax or another function)')
        exit(1)

        
def lazy_integral(f,a,b,tolerance,nmax):
    """
    Functions that integrate f from the intevral [a,b] using the dynamic simpson's rule above. 
    It basically sums over the value at each subinterval and returns the result and the number
    of counts it took for the integral to converge
        -Inputs: f (function) Function that we are interested in integrating 
                 a (float) Lower bound of the integral
                 b (float) Upper bound of the integral
                 tolerance (float) Tolerance on the error of the integral 
                 nmax (int) Maximum number of steps for the recursion. If reached, we kill
                                       it
        -Outputs: I (float) integral value 
                  count (int) Number of steps for the recursion to converge
    """

    fa = f(a)
    fb = f(b)
    guess, h, fh = simpson(f,a,b,fa,fb)
    results = np.array(dynamic_simpson(f,a,b,h,fa,fb,fh,guess,tolerance,nmax))
    count    = np.sum(results[1::2])
    integral = np.sum(results[0::2])
    return integral, count

def integrate(functions, x_min,x_max,tolerance = 1e-7,nmax=1000):
    """
    Integrates different inputed functions
        -Inputs: funcs (array of functions). These are the functions to be integrated
                 x_min (float) Lower bound of the integral 
                 x_max (float) Upper bound of the integral
                 tolerance (float) Set to default to 1e-7, this is the tolerance of the integral
                 nmax (float) Maxium number of steps given fo the integral to diverge. 
    """
    for func in functions: 
        I, count = lazy_integral(func[0],x_min,x_max,tolerance,nmax)
        error = np.abs(I - func[2](x_max)+func[2](x_min))
        return I, error