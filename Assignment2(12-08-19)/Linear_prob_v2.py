import math
import time

a = [[2,1,1,3,2],
[1,2,2,1,1],
[1,2,9,1,5],
[3,1,1,7,1],
[2,1,5,1,8]]

b = [-2,4,3,-5,1]

def disp_mat(z,n):
    for row in range(n):
        print(z[row])

def doolittle(x,b):
    for i in range(0,len(x)):
        if (len(x[i])==len(x)):
            pass
        else:
            print("\n Non-Square matrix, returning None")
            return(None)

    #for square matrix
    n = len(x)
    u= [[0 for i in range(n)] for j in range(n)]
    l= [[0 for i in range(n)] for j in range(n)]

    for i in range(n):
        l[i][i]=1
        for j in range(i-1):
            s= sum(l[i][k]*u[k][j] for k in range(j-1))
            l[i][j] = (x[i][j] - s) / u[j][j]
        for j in range(i,n):
            s= sum(u[k][j] * l[i][k] for k in range(j-1))
            u[i][j] = x[i][j] - s
    print("\n doolittle U")
    disp_mat(u,n)
    print("\n doolittle L")
    disp_mat(l,n)

    # two steps : 1) LZ=B  2)UX=Z

    z=[0 for i in range(n)]
    sol=[0 for i in range(n)]

    z=[0]*n
    sol = [0]*n

    for i in range(n):
        s= sum(l[i][j]*z[j] for j in range(i-1))
        z[i] = b[i]- s

    for c in range(n):
        i=(n-1)-c
        s=sum(u[i][j]*sol[j] for j in range(i+1,n))
        sol[i]= (z[i]-s)/u[i][i]

    print(sol)
    return(sol)


def croute(x,b):
    for i in range(0,len(x)):
        if (len(x[i])==len(x)):
            pass
        else:
            print("\n Non-Square matrix, returning None")
            return(None)

    #for square matrix
    n = len(x)
    u= [[0 for i in range(n)] for j in range(n)]
    l= [[0 for i in range(n)] for j in range(n)]

    for i in range(n):
        u[i][i]=1
        for j in range(n):
            s= sum(l[j][k]*u[k][i] for k in range(i-1))
            l[j][i]= (a[j][i]-s)

        for j in range(i+1,n):
            s= sum(l[i][k]*u[k][j] for k in range(i-1))
            u[i][j] = (a[i][j] - s)/l[i][i]

    print("\n croute U")
    disp_mat(u,n)
    print("\n croute L")
    disp_mat(l,n)

    z=[0 for i in range(n)]
    sol=[0 for i in range(n)]

    for i in range(n):
        s= sum(l[i][j]*z[j] for j in range(i-1))
        z[i] = (b[i]- s)/l[i][i]

    for c in range(n):
        i=(n-1)-c
        s=sum(u[i][j]*sol[j] for j in range(i+1,n))
        sol[i]= (z[i]-s)

    print(sol)
    return(sol)

def cholesky(x,b):
    n = len(x)
    for i in range(0,len(x)):
        if (len(x[i])==len(x)):
            pass
        else:
            print("\n Non-Square matrix, returning None")
            return(None)
    for i in range(n):
        for j in range(i):
            if(x[i][j] != x[j][i]):
                print("\n Non-Symmetric matrix, returning None")
                return(None)

    u= [[0 for i in range(n)] for j in range(n)]
    z=[0 for i in range(n)]
    sol=[0 for i in range(n)]

    for i in range(n):
        s = sum(u[i][k]*u[i][k] for k in range(i-1))
        u[i][i] = math.sqrt(x[i][i]-s)

        for j in range(i+1,n):
            s = sum(u[k][i]*u[k][j] for k in range(i-1))
            u[i][j] = (x[i][j]-s)/u[i][i]

    print("\n cholesky U")
    disp_mat(u,n)

    for i in range(n):
        s = sum(u[i][j]*z[j] for j in range(i-1))
        z[i]=(b[i]-s)/u[i][i]

    for c in range(n):
        i=(n-1)-c
        s= sum(u[i][j]*sol[j] for j in range(i+1,n))
        sol[i]= (z[i]-s)/u[i][i]

    print(sol)
    return(sol)

ans= cholesky(a,b)
ans = doolittle(a,b)
ans = croute(a,b)

# try:
#     ans= cholesky(a,b)
#     if(ans!= None):
#         pass
# except:
#     try:
#         ans = doolittle(a,b)
#     except:
#         ans = croute(a,b)


print ("\n CPU time: ", time.process_time(),'s')
