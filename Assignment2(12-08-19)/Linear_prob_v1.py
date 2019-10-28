
a = [[2,1,1,3,2],
[1,2,2,1,1],
[1,2,9,1,5],
[3,1,1,7,1],
[2,1,5,1,8]]

b = [-2,4,3,-5,1]

# by Doolittle method: L has 1 as its diagonal elements

def disp_mat(z,n):
    for row in range(n):
        print(z[row])

def doolittle(x):
    for i in range(0,len(x)):
        if (len(x[i])==len(x)):
            pass
        else:
            print("\n Non-Square matrix, returning (0,0)")
            return(0,0)

    #for square matrix
    n = len(x)
    u= [[0 for i in range(n)] for j in range(n)]
    l= [[0 for i in range(n)] for j in range(n)]
    t= [[0 for i in range(n)] for j in range(n)]
    m= [[0 for i in range(n)] for j in range(n)]

    u[0][0]=1
    print(u)

    for i in range(n):
        l[i][i]=1
        for j in range(i-1):
            s=0
            for k in range(j-1):
                s = s+ u[k][j] * l[i][k]
            l[i][j] = (x[j][i] - s) / u[j][j]
            print("L[",i,"][",j,"] :",l[i][j])
            m[i][j]= l[i][j]
            
        for j in range(i,n):
            s=0
            for k in range(j-1):
                s = s + u[k][j] * l[i][k] 
            u[i][j] = x[i][j] - s
            print("U[",i,"][",j,"] :",u[i][j])
            t[i][j]= u[i][j]
        
    #print(m)
    #print(p)
    disp_mat(l,n)
    disp_mat(u,n)
    return(l,u)

def croute(x):
    for i in range(0,len(x)):
        if (len(x[i])==len(x)):
            pass
        else:
            print("\n Non-Square matrix, returning (0,0)")
            return(0,0)

    #for square matrix
    n = len(x)
    u= [[0]*n]*n
    l= [[0]*n]*n

    for i in range(n):
        u[i][i]=1
        for j in range(i-1):
            s1= sum(l[j][k]*u[k][i] for k in range(i-1))
            l[j][i]= (a[j][i]-s1)
            print("L[",i+1,"][",j+1,"] :",l[i][j])
        for j in range(i+1,n):
            s1= sum(l[i][k]*u[k][j] for k in range(1,i-1))
            u[i][j] = (a[i][j] - s1)/l[i][i]
            print("U[",i+1,"][",j+1,"] :",u[i][j])
            
    disp_mat(l,n)
    disp_mat(u,n)
    return(l,u)

#def cholesky(x):
    


            
#def LU_solve(l,u,b):
    


(l,u)=doolittle(a)
#(l,u)=croute(a)
