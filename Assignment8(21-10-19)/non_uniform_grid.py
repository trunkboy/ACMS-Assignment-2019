import time
import matplotlib.pyplot as plt

def linearsolver(A,b):
  n = len(A)
  M = A

  i = 0
  for x in M:
   x.append(b[i])
   i += 1

  for k in range(n):
   for i in range(k,n):
     if abs(M[i][k]) > abs(M[k][k]):
        M[k], M[i] = M[i],M[k]
     else:
        pass

   for j in range(k+1,n):
       q = float(M[j][k]) / M[k][k]
       for m in range(k, n+1):
          M[j][m] -=  q * M[k][m]

  x = [0 for i in range(n)]

  x[n-1] =float(M[n-1][n])/M[n-1][n-1]
  for i in range (n-1,-1,-1):
    z = 0
    for j in range(i+1,n):
        z = z  + float(M[i][j])*x[j]
    x[i] = float(M[i][n] - z)/M[i][i]
  return x

def mat_mul(a,b,x):
    n = len(x)

    sol = [0 for i in range(n)]
    err = [0 for i in range(n)]

    for i in range(n):
        sum=0
        for j in range(n):
            sum = sum+a[i][j]*x[j]
        sol[i] = sum

    for i in range(n):
        err[i] = b[i]-sol[i]

    return(err)

def thomas(a,b):
    n= len(b)
    for i in range(1,n):
        a[i][i-1]=a[i][i-1]/a[i-1][i-1]
        a[i][i] = a[i][i]-a[i][i-1]*a[i-1][i];
        b[i]= b[i]-a[i][i-1]*b[i-1]
        a[i][i-1]=0

    #backward substitution
    x=[0 for i in range(n)]
    x[n-1]= b[n-1]/a[n-1][n-1]
    for k in range(n-1):
        i=n-2-k
        x[i] = ( b[i] + x[i+1] )/ a[i][i]

    #print(x)
    return(x)

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

def grid_1d_gp(r,n):
    x = [0 for i in range(n+1)]
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

def grid_1d_gp2(r,n):
    x = [0 for i in range(n+1)]
    rng=(0,1)

    #gp sum
    sum=0
    for i in range(n):
        sum = sum + pow(r,i)
    x0 = (rng[1]-rng[0])/sum
    x[0] = rng[0]
    x[n] = rng[1]

    step = [ 0 for i in range(n)]
    for i in range(n):
        step[i] = x0*pow(r,i)

    for i in range(n-1,0,-1):
        x[i] = x[i+1] - step[n-1-i]
    print(x)
    return(x)

def matrix_form(x,bc):

    p=1
    u=1
    t=0.02

    n = len(x)-2
    a = [[0 for i in range(n)] for j in range(n)]

    for j in range(2,n):
        i=j-1

        beta = ( (-p*u)/(x[i+1]-x[i-1]) ) - ( (2*t)/((x[i]-x[i-1])*(x[i+1]-x[i-1])) )
        alpha = ( (2*t)/((x[i+1]-x[i-1])*(x[i+1]-x[i])) ) - ( (2*t)/((x[i+1]-x[i-1])*(x[i]-x[i-1])) )
        gamma = ( (p*u)/(x[i+1]-x[i-1]) ) - ( (2*t)/((x[i+1]-x[i])*(x[i+1]-x[i-1])) )

        a[i][i-1] = gamma
        a[i][i] = alpha
        a[i][i+1] = beta

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

def solve_plot( r,n):
    bc= (0,1)
    x = grid_1d_gp(r,n)

    a,b = matrix_form(x,bc)

    phi = gauss_elm(a,b)
    # phi = thomas(a,b)
    # phi = linearsolver(a,b)

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
