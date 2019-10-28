import numpy as np
import matplotlib.pyplot as plt
import time

r,n = 0.3,50

p=1
u=1
t=0.02

def grid(r,n):
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

x = grid(r,n)
a = np.zeros(n+1)
b = np.zeros(n+1)
c = np.zeros(n+1)
d = np.zeros(n+1)
phi = np.zeros(n+1)

#boundary conditions
bc = (0,1)
phi[0] = bc[0]
phi[n] = bc[1]

#Calculating coefficients
for i in range(1,n):
    gamma = ( (-p*u)/(x[i+1]-x[i-1]) ) - ( (2*t)/((x[i]-x[i-1])*(x[i+1]-x[i-1])) )
    alpha = ( (2*t)/( (x[i+1]-x[i-1])*(x[i+1]-x[i]) ) ) + ( (2*t)/((x[i+1]-x[i-1])*(x[i]-x[i-1])) )
    beta = ( (p*u)/(x[i+1]-x[i-1]) ) - ( (2*t)/((x[i+1]-x[i])*(x[i+1]-x[i-1])) )

    if(i is 1):
        b[i]=alpha
        c[i] =beta
        d[i] = 0 - (gamma*bc[0])
    elif(i is n-1):
        b[i]=alpha
        a[i]=gamma
        d[i] =  0 - (beta*bc[1])
    else:
        a[i] = gamma
        b[i] = alpha
        c[i] = beta

#Thomas algorithm
for i in range(2,n):
    a[i] = a[i]/b[i-1]
    b[i] = b[i] - a[i]*c[i-1]
    d[i] = d[i] - a[i]*d[i-1]

phi[n-1] = d[n-1] / b[n-1]
for i in range(n-2,0,-1):
    phi[i]= (d[i] - c[i]*phi[i+1])/b[i]

print(phi)
print(x)
plt.plot(x, phi,'r+')
plt.title("phi vs. x")
plt.ylabel('phi')
plt.xlabel('x')
plt.show()
