import time
import matplotlib.pyplot as plt

def gauss_siedel(a,b):
    def diagonal_dominance(a):
        n=len(a)
        dom=0
        for i in range(n):
            dom=0
            for j in range(n):
                if(i!=j and abs(a[i][i])>abs(a[i][j])):
                    dom=dom+1;
            if(dom == (n-1)):
                print("\ndiagionally dominant")
                return(1)
            else:
                print("checking for next row")

        print("Not diagionally dominant, needs reordering")
        return(0)

    #check diagonal dominance:
    if(not(diagonal_dominance(a))):
        print("NA_returning 0 ")
        return(0)

    n=len(a)
    guess = [0 for i in range(n)]
    lim_err = 0.0001
    #forming equation matrix
    z=[[a[i][j] for i in range(n)] for j in range(n)]
    for i in range(n):
       z[i][i] = 0

    x=guess[:]
    x_prev = x[:]
    err = [0]*n
    iteration = 0
    while(1):
        #calculating iteration with error
        max_err=0
        x_prev = x[:]
        for i in range(n):
            s = sum(z[i][j]*x[j] for j in range(n))
            x[i]=(b[i]-s)/a[i][i]

            #error
            err[i] = ((x[i]-x_prev[i])/(x[i]))
            if(err[i]<0):
                err[i] = (-1)*err[i]
            if(max_err<err[i]):
                max_err=err[i]

        #print("\n\nError is :",err)
        #print("x_prev :",x_prev)
        #print("x :",x)
        iteration = iteration+1
        if(max_err<lim_err):
            break
    # print("\n\nIterations:",iteration)
    # print("solution:",x)
    # print("initial guess was:",guess)
    # print ("\n CPU time: ", time.process_time(),'s')
    return(x)


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

# formulate tridiagonal matrix

def solve_plot(delx,scheme):
    def matrix_form(delx,scheme):
        rng = 1.0
        n = int((rng/delx) -1)
        #print (n)
        a = [[0 for i in range(n)] for j in range(n)]
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
            a[i][i-1] = gamma
            a[i][i] = beta
            a[i][i+1] = alpha
        #row 0
        a[0][0]=beta
        a[0][1] =alpha

        #row n-1
        a[n-1][n-1]=beta
        a[n-1][n-2]=gamma

        #rhs vector
        b = [0 for i in range(n)]
        b[0] = -gamma
        b[n-1] =  0
        return(a,b,xm,n)

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

    #display matrices and check
    a,b,xm,n = matrix_form(delx,scheme)
    # for i in range(n):
    #     print(a[i])
    # for i in range(n):
    #     print(b[i])

    # phi=thomas(a,b)
    # phi=gauss_elm(a,b)
    # phi=linearsolver(a,b)
    phi=gauss_siedel(a,b)

    mod_phi = [0 for i in range(n+1)]
    for i in range(len(phi)):
        mod_phi[i+1] = phi[i]
    mod_phi[0] = 1
    mod_phi[n]= 0
    print(mod_phi)
    print(xm)
    print ("\n CPU time: ", time.process_time(),'s')

    #plotting
    plt.plot(xm,mod_phi)
    plt.xlabel('x')
    plt.ylabel('phi')
    meta = str(delx)
    meta = "phi vs. x : dx = " + meta
    plt.title(meta)
    plt.show()

# solve_plot(0.1,"cds")

solve_plot(0.01,"cds")

# solve_plot(0.01,"fds")

# solve_plot(0.01,"bds")
