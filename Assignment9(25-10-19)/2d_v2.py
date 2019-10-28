=import time
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np


def uniform_grid(lims, n):
    rng = lims[1]-lims[0]
    delx=rng/n

    xm = [0 for i in range(n+1)]
    for i in range(n):
        xm[i+1] = xm[i]+delx

    return(xm,delx)

def matrix_form_solve_plot(lx,nx,ly,ny):
    global x,dx,y,dy
    x,dx = uniform_grid(lx,nx)
    y,dy = uniform_grid(ly,ny)

    e= 1
    c= -2*(1-(dx*dx)/(dy/dy))
    w= 1

    s= (dx*dx)/(dy*dy)
    n= (dx*dx)/(dy*dy)

    lx = len(x)-1
    ly = len(y)-1

    global t
    t = [ [0 for i in range(lx+1)] for j in range(ly+1) ]
    nt = [ [0 for i in range(lx+1)] for j in range(ly+1) ]
    Qi =0

    # boundary conditions
    # 1. left
    for i in range(ly+1):
        t[i][0] = 75
    # 2. Right
    for i in range(ly+1):
        t[i][lx] = 100
    # 3. Top
    for i in range(lx+1):
        t[ly][i] = 300
    # 4. Bottoms
    for i in range(lx+1):
        t[0][i] = 50


    #jacobi run
    run = 0
    while(run<50):
        for i in range(1,ly):
            for j in range(1,lx):
                nt[i][j] = (Qi - ( w*t[i-1][j] + e*t[i+1][j] + s*t[i][j]+ n*t[i][j] ))/(c)

        #reiterate
        run = run+1
        for i in range(1,ly):
            for j in range(1,lx):
                t[i][j]= nt[i][j]

    #plotting contour
    X, Y = np.meshgrid(x, y)
    fig,ax=plt.subplots(1,1)
    cp = ax.contourf(X, Y, t)
    fig.colorbar(cp) # Add a colorbar to a plot
    ax.set_title('Filled Contours Plot')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    plt.show()

    #plotting 3D
    X, Y = np.meshgrid(x, y)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.contour3D(X, Y, t, 50, cmap='binary')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('t')
    ax.set_title('3D contour')
    plt.show()

nx=10
ny=10

limx = (0,2.4)
limy = (0,3.0)

matrix_form_solve_plot(limx,nx,limy,ny)
