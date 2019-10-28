import time
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np

# A roll of steel 100C at left end and 25C at right end, length = 0.05m
# find temperature variation from t=0 to t=9s

def uniform_grid(lims, n):
    rng = lims[1]-lims[0]
    delx=rng/n

    xm = [0 for i in range(n+1)]
    for i in range(n):
        xm[i+1] = xm[i]+delx

    return(xm,delx)

def disp(x):
    for i in range(len(x)):
        print(x[i])

def gauss_elm(A,B):
    n= len(B)
    # step 1: Gaussian elimination.
    i=0
    while i < n:
        # pivots
        pivot = A[i][i]
        j=i+1
        while j<n:
            r = A[j][i]/pivot
            # row opreation
            k=i
            while k<n:
                A[j][k] = A[j][k] - A[i][k]*r
                k=k+1

            B[j]=B[j]-B[i]*r
            j=j+1
        i=i+1

    #Back Substitution from nth row
    x= [0 for i in range(n)]

    i = n-1
    x[i] = B[i]/A[i][i]
    i=i-1
    while i>=0:
        sum = 0
        k=i+1
        while k<n:
            sum = sum + A[i][k]*x[k]
            k=k+1
        x[i]=(B[i]-sum)/A[i][i]
        i=i-1
    return(x)

def plot_temp(x,time,temp):
    #plotting filled contour
    X, Y = np.meshgrid(x, time)
    fig,ax=plt.subplots(1,1)
    cp = ax.contourf(X, Y, temp,100,cmap = 'viridis')
    fig.colorbar(cp) # Add a colorbar to a plot
    ax.set_title('Filled Contour Plot')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('time (s)')
    plt.show()

def temperature_array_explicit(nx):
    bcx = (100,25)
    bct = (0,9)
    alpha = 540/(7800*490)

    x,dx = uniform_grid(( 0,0.05 ), nx)

    dt=(dx*dx)/(4*2*alpha)
    nt = int((bct[1] - bct[0])/dt)
    print("nt = ",nt)

    time,dt = uniform_grid(( 0,9 ), nt)
    temp = [[20 for i in range(nx+1)] for j in range(nt+1)]

    #boundary conditions
    for i in range(nt+1):
        temp[i][0] = bcx[0]
        temp[i][nx] = bcx[1]

    #ntemp = [[temp[i][j] for i in range(nx+1)] for j in range(nt+1)]
    r=alpha*dt/(dx*dx)


    #EXPLICIT SOLVING
    # dt<(dx^2)/(4*alpha)
    #constants
    e = r
    w = e
    p = 1-2*r

    #Iterations
    for j in range(nt):
        for i in range(1,nx):
            temp[j+1][i] = e*temp[j][i+1] + p*temp[j][i] + w*temp[j][i-1]

    plot_temp(x,time,temp)

def temperature_array_implicit(nx,nt):
    bcx = (100,25)
    bct = (0,9)
    alpha = 540/(7800*490)

    x,dx = uniform_grid(( 0,0.05 ), nx)
    time,dt = uniform_grid(( 0,9 ), nt)
    temp = [[20 for i in range(nx+1)] for j in range(nt+1)]

    #boundary conditions
    for i in range(nt+1):
        temp[i][0] = bcx[0]
        temp[i][nx] = bcx[1]

    #ntemp = [[temp[i][j] for i in range(nx+1)] for j in range(nt+1)]
    r=alpha*dt/(dx*dx)

    #IMPLICIT SOLVING
    e = -r
    p = 2*r+1
    w = e

    for k in range(nt):
        a = [[0 for i in range(nx-1)] for j in range(nx-1)]
        b = [0 for i in range(nx-1)]
        for i in range(nx-1):
            Qi = temp[k][i+1]
            if(i is 0):
                a[0][0]=p
                a[0][1] =e
                b[0] = Qi - (w*bcx[0])
            elif(i is nx-2):
                a[i][i]=p
                a[i][i-1]=w
                b[i] = Qi - (e*bcx[1])
            else:
                a[i][i-1] = w
                a[i][i] = p
                a[i][i+1] = e
                b[i] = Qi

        m = [[a[i][j] for i in range(nx-1)] for j in range(nx-1)]
        n = [b[i] for i in range(nx-1)]

        zeta = gauss_elm(m,n)
        for i in range(0,nx-1):
            temp[k+1][i+1] = zeta[i]
    plot_temp(x,time,temp)

def temperature_array_cr(nx):
    bcx = (100,25)
    bct = (0,9)
    alpha = 540/(7800*490)

    x,dx = uniform_grid(( 0,0.05 ), nx)

    dt=(dx*dx*0.8)/(4*alpha)
    nt = int((bct[1] - bct[0])/dt)
    print("nt = ",nt)

    time,dt = uniform_grid(( 0,9 ), nt)
    temp = [[20 for i in range(nx+1)] for j in range(nt+1)]

    #boundary conditions
    for i in range(nt+1):
        temp[i][0] = bcx[0]
        temp[i][nx] = bcx[1]

    #ntemp = [[temp[i][j] for i in range(nx+1)] for j in range(nt+1)]
    r=alpha*dt/(dx*dx)
    #CRANK -NICHOLSON
    e = -r/2
    p = r+1
    w = e

    for k in range(nt):
        a = [[0 for i in range(nx-1)] for j in range(nx-1)]
        b = [0 for i in range(nx-1)]
        for i in range(nx-1):
            Qi = (1-r)*temp[k][i+1] + r*temp[k][i+2]/2  + r*temp[k][i]/2

            if(i is 0):
                a[0][0]=p
                a[0][1] =e
                b[0] = Qi - (w*bcx[0])
            elif(i is nx-2):
                a[i][i]=p
                a[i][i-1]=w
                b[i] = Qi - (e*bcx[1])
            else:
                a[i][i-1] = w
                a[i][i] = p
                a[i][i+1] = e
                b[i] = Qi

        m = [[a[i][j] for i in range(nx-1)] for j in range(nx-1)]
        n = [b[i] for i in range(nx-1)]

        zeta = gauss_elm(m,n)
        for i in range(0,nx-1):
            temp[k+1][i+1] = zeta[i]
    plot_temp(x,time,temp)

nx=100
nt=100

# temperature_array_explicit(nx)
# temperature_array_implicit(nx,nt)
temperature_array_cr(20)
