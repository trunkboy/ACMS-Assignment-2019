import time

# a= [[2.08,-1,0,0],
# [-1,2.08,-1,0],
# [0,-1,2.08,-1],
# [0,0,-1,2.08]]


# b= [41.6,1.6,1.6,201.6]


a = [2.08,2.08,2.08,2.08]
#super diagonal
b = [-1,-1,-1]
#sub diagonal
c = [-1,-1,-1]
#result vector
d= [41.6,1.6,1.6,201.6]
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

    print(x)

thomas(a,b,c,d)
print ("\n CPU time: ", time.process_time(),'s')
