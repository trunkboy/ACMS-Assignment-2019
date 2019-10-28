import time
import matplotlib.pyplot as plt
import numpy as np

def thomas2(a,b,c,d):
    n = len(a)-1
    for i in range(2,n):
        a[i] = a[i]/b[i-1]
        b[i] = b[i] - a[i]*c[i-1]
        d[i] = d[i] - a[i]*d[i-1]

    phi[n-1] = d[n-1] / b[n-1]
    for i in range(n-2,0,-1):
        phi[i]= (d[i] - c[i]*phi[i+1])/b[i]

    return(phi)

def grid_1d_gp(r,n):
    x = np.zeros(n+1)
    rng=(0,1)

    #gp sum
    sum=0
    for i in range(n):
        sum = sum + pow(r,i)
    x0 = (rng[1]-rng[0])/sum
    x[0] = rng[0]

    for i in range(n):
        x[i+1] = x[i] + x0*pow(r,i)
    # print(x)
    return(x)

def matrix_form(x,bc):

    p=1
    u=1
    t=0.02

    n = len(x)-1
    print(n)
    a = np.zeros(n+1)
    b = np.zeros(n+1)
    c = np.zeros(n+1)
    d = np.zeros(n+1)

    for i in range(2,n-1):
        beta = ( (-p*u)/(x[i+1]-x[i-1]) ) - ( (2*t)/((x[i]-x[i-1])*(x[i+1]-x[i-1])) )
        alpha = ( (2*t)/((x[i+1]-x[i-1])*(x[i+1]-x[i])) ) - ( (2*t)/((x[i+1]-x[i-1])*(x[i]-x[i-1])) )
        gamma = ( (p*u)/(x[i+1]-x[i-1]) ) - ( (2*t)/((x[i+1]-x[i])*(x[i+1]-x[i-1])) )

        a[i] = gamma
        b[i] = alpha
        c[i] = beta

    #row 0
    a[0][0]=alpha
    a[0][1] =beta

    #row n-1
    a[n-1][n-1]=alpha
    a[n-1][n-2]=gamma

    #rhs vector
    b = [0 for i in range(n)]
    b[0] = 0 - (gamma*bc[0])
    b[n-1] =  0 - (beta*bc[1])

    return(a,b)

def solve_plot(r,n):
    bc= (0,1)
    x = grid_1d_gp(r,n)

    a,b = matrix_form(x,bc)

    phi = gauss_elm(a,b)
    # phi = thomas(a,b)
    # phi = linearsolver(a,b)
    # phi = thomas(a,b)

    print('err is:')
    print(mat_mul(a,b,phi))

    mod_phi = [0 for i in range(n+1)]
    for i in range(len(phi)):
        mod_phi[i+1] = phi[i]
    mod_phi[0] = bc[0]
    mod_phi[n]= bc[1]
    print('phi is:')
    print(mod_phi)
    print('x is:')
    print(x)
    print ("\n CPU time: ", time.process_time(),'s')

    #plotting
    plt.plot(x,mod_phi,'r+')
    plt.xlabel('x')
    plt.ylabel('phi')
    plt.title('phi vs. x')
    plt.show()

#main

solve_plot(0.70,50)
