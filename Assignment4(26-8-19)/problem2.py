# THe manning equation can be written for a rectangular open channel flow as
# Q = sqrt(S) * (BH)^(5/3)/(n*(B+2H)^(2/3))
# where Q is flow rate, S is slope and H is depth, n is manning rougness coefficient, b is breadth
# Develop fixed point iteration snippet to solve for H, given:
# Q=5 ; s = 0.0002, B=20m ,n=0.03 ; wiht error limit of 0.05%
# Prove that your scheme coverges for all initial guess greater than or equal zero.

import time
import math

Q=5
s=0.0002
B=20
n=0.03

init =25
lim = 0.0005

def manning(Q,B,n,s,H):
    a=pow((B+2*H),(2/3))
    z = pow(s,0.5)
    f = (Q*n*a)/z
    ret = pow(f,(3/5))*(1/B)
    return(ret)

# def manning_H2(Q,B,n,s,H):
#     a = pow(B*H,(5/3))*sqrt(s)
#     z = pow((a/(n*Q)),(3/2))
#     ret = 0.5*(B-z)
#     return(ret)

#call these two possible functions simultaneously to find the root of the function

def fpi(f,args,init_guess,err_lim):
    #using function manning
    itr = 0
    x = init_guess
    print("\ninitial value :",x)
    while(1):
        itr = itr+1
        argm = args + (x,)
        y = f(*argm)
        err = (y-x)/y
        if(err<0):
            err = (-1)*err
        print("\nx :",y)
        print("error :",err)
        if(err<err_lim):
            break
        x=y
    print("\n Iterations:",itr)
    print ("\n CPU time: ", time.process_time(),'s')
    return(y)

ans = fpi(manning,(Q,B,n,s),init,lim)
