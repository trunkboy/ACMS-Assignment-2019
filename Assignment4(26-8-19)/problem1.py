import time
from math import pi as pi

#  problem statement: To find the maximum of a function given as:

# A spherical tank is to be designed to hold water for a small village in a developing country.
# The volume of liquid it can hold can be computed as :  V=(pi*h^2 * (3R-h)/3)
# with R=3, what depth must the tank be filled so that it holds 30m^3.
# Solve the problem by developing a code using NR method and compute your answer for 5 decimal points accuracy.
# Calc the number of iterations required to converge to the given criteria.
# Compute approx relative error after every iteration.

err_lim = 0.00001
init_guess= 3
r_val=3

def V(h,r):
    return(pi*h*h*(3*r-h)/3)

def dV(h,r):
    return(pi*h*(2*r-h))

def NR(r,vol,guess_h,err_lim):
    x= guess_h
    err = 1
    itr = 0
    
    while(err>err_lim):
        itr = itr+1
        y= x-( (V(x,r)-vol)/(dV(x,r)) )
        err = (y-x)/y
        if(err<0):
            err = err*(-1)
        print("\n X[",itr-1,"] =",x,"\n X[",itr,"] =",y)
        print("Error is :",err)
        x=y

    print("\n Solution :",y)
    print("Iterations :",itr)
    print ("\n CPU time: ", time.process_time(),'s')
    return(y,itr)

z = NR(r_val,30,init_guess,err_lim)
