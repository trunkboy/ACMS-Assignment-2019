import time
import matplotlib.pyplot as plt
import numpy as np
from numba import jit, f8

## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
@jit(f8[:] (f8[:],f8[:],f8[:],f8[:] ))
def TDMAsolver(b, c, a, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    and to http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
    '''
    nf = len(d) # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
    for it in range(1, nf):
        mc = ac[it-1]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1]
        dc[it] = dc[it] - mc*dc[it-1]

    xc = bc
    xc[-1] = dc[-1]/bc[-1]

    for il in range(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

    return xc

def thomas(a,b,c,d):
    n= len(d)
    for i in range(1,n):
        c[i-1]=c[i-1]/a[i-1]
        a[i] = a[i]-c[i-1]*b[i-1];
        d[i]= d[i]-c[i-1]*d[i-1]
        c[i-1]=0

    #backward substitution
    x=[0 for i in range(n)]
    x[n-1]= d[n-1]/a[n-1]
    for k in range(n-1):
        i=n-2-k
        x[i] = ( d[i] + x[i+1] )/ a[i]

    #print(x)
    return(x)

def TDMASolve(b, c, a, d):
    n = len(a)
    ac, bc, cc, dc = map(np.array, (a, b, c, d))
    xc = []
    for j in range(2, n):
        if(bc[j - 1] == 0):
            ier = 1
            return
        ac[j] = ac[j]/bc[j-1]
        bc[j] = bc[j] - ac[j]*cc[j-1]
    if(b[n-1] == 0):
        ier = 1
        return
    for j in range(2, n):
        dc[j] = dc[j] - ac[j]*dc[j-1]
    dc[n-1] = dc[n-1]/bc[n-1]
    for j in range(n-2, -1, -1):
        dc[j] = (dc[j] - cc[j]*dc[j+1])/bc[j]
    return dc

def solve_plot(delx,scheme):

    def matrix_form(delx,scheme):
        rng = 1.0
        n = int((rng/delx) -1)
        #print (n)
        a = [0 for i in range(n)]
        b = [0 for i in range(n-1)]
        c = [0 for i in range(n-1)]

        t = 0.1

        xm = [0 for i in range(n+1)]
        for i in range(n):
            xm[i+1] = xm[i]+delx
        #print (x)
        if(scheme is "cds"):
            alpha = 1+ (delx/2*t)
            beta = -2-(delx*delx/t)
            gamma = 1- (delx/2*t)
        elif(scheme is "bds"):
            alpha = 1
            beta = -2-(delx*delx/t) +(delx/t)
            gamma = 1- (delx/t)
        elif(scheme is "fds"):
            alpha = 1+ (delx/t)
            beta = -2-(delx*delx/t) -(delx/t)
            gamma = 1

        for i in range(1,n-1):
            c[i-1] = gamma
            a[i] = beta
            b[i] = alpha

        #row 0
        a[0]=beta
        b[0] = alpha

        #row n-1
        a[n-1]=beta
        c[n-2]=gamma

        #rhs vector
        d = [0 for i in range(n)]
        d[0] = -gamma
        d[n-1] =  0

        # print(a)
        # print(b)
        # print(c)
        # print(d)
        return(a,b,c,d,xm,n)

    #display matrices and check
    a,b,c,d,xm,n = matrix_form(delx,scheme)
    # for i in range(n):
    #     print(a[i])
    # for i in range(n):
    #     print(d[i])


    #phi=thomas(a,b,c,d)
    # phi=TDMAsolver(a,b,c,d)
    phi = TDMASolve(a,b,c,d)


    mod_phi = [0 for i in range(n+1)]
    for i in range(len(phi)):
        mod_phi[i+1] = phi[i]
    mod_phi[0] = 1
    mod_phi[n]= 0
    print(mod_phi)
    print(xm)
    print ("\n CPU time: ", time.process_time(),'s')

    #plotting
    plt.plot(mod_phi, xm)
    plt.xlabel('x')
    plt.ylabel('phi')
    meta = str(delx)
    meta = "phi vs. x : dx = " + meta
    plt.title(meta)
    plt.show()


# solve_plot(0.1,"cds")


solve_plot(0.025,"cds")

# solve_plot(0.01,"fds")

# solve_plot(0.01,"bds")
