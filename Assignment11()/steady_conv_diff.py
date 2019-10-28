import time
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

def disp(x):
    for i in range(len(x)):
        print(x[i])

def matrix_form(nx,ny):
    limx = (0,1)
    limy = (0,1)

    x,dx = uniform_grid(limx,nx)
    y,dy = uniform_grid(limy,ny)

    # convection velocity matrix_form
    u = [ [0 for i in range(nx+1)] for j in range(ny+1) ]
    v = [ [0 for i in range(nx+1)] for j in range(ny+1) ]

    for j in range(0,ny+1):
        for i in range(0,nx+1):
            #west
            u[j][i] = x[i]
            #east
            try:
                u[j][i+1] = x[i+1] #IS this required?
            except:
                pass
            #south
            v[j][i] = -y[j]


    #diffusion and convection terms coefficients
    for j in range(0,ny+1):
        for i in range(0,nx+1):
            
