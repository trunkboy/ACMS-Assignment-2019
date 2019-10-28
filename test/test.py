import time
import math

#globals


g = 9.81
Q = 20

def fun(g_,Q_,y):
    a = g_*pow(y,3)
    b = pow ((6+y),3)
    c = (3+y)*Q_*Q_

    fun_val = ((a*b)/8)-c

    return(fun_val)

def bisect_root(guess,g,Q):
    a= guess[0]
    b= guess[1]
    itr = 0
    while(1):
        f_a = fun(g,Q,a)
        f_b = fun(g,Q,b)

        m = (a+b)/2
        f_m = fun(g,Q,m)

        if(f_a*f_m > 0):     #same signs so they should replace
            a = m
        elif(f_b*f_m > 0 ):
            b = m
        else:
            print("\n Problem with elif")

        if(itr):
            #erroe
            err = (m - m_prev)/m
            if(err<0):
                err=-1*err
            print("\n Error is :",err)
            if(err<0.001):
                print("\n Process over")
                break
        m_prev=m
        itr = itr +1

    print("Solution is:", m )
    print("Error is:",err)
    print("function value is at the solution is:",fun(g,Q,m))
    print("Number of iterations",itr)

bisect_root((0.5,2.5),g,Q)
