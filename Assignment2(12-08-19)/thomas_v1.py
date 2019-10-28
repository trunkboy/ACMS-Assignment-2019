import time

a= [[2.08,-1,0,0],
[-1,2.08,-1,0],
[0,-1,2.08,-1],
[0,0,-1,2.08]]

b= [41.6,1.6,1.6,201.6]

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

    print(x)
    return(x)

z=thomas(a,b)
print('err is:',mat_mul(a,b,z))
print ("\n CPU time: ", time.process_time(),'s')
