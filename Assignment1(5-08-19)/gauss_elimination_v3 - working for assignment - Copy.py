import time

#For matrix system AX=B, solve using gaussion elimination and produce the output along with computational time

A = [[25,5,1],
    [64,8,1],
    [144,12,1]]

B = [106.8,177.2,279.2]

# Step 1: Gaussian elimination.
i=0
while i < 3:
    # pivots
    pivot = A[i][i]
    j=i+1
    while j<3:
        r = A[j][i]/pivot
        # row opreation
        k=i
        while k<3:
            A[j][k] = A[j][k] - A[i][k]*r
            k=k+1

        B[j]=B[j]-B[i]*r
        #print (A)
        #print (B)
        #print("\n")

        j=j+1

    i=i+1

#print (A)
#print (B)

#Step 2: Back Substitution from nth row

x= [0]*3
n=3
i = n-1
x[i] = B[i]/A[i][i]
i=i-1
while i>=0:
    sum = 0
    k=i+1
    while k<n:
        #print("kth")
        sum = sum + A[i][k]*x[k]
        k=k+1
    x[i]=(B[i]-sum)/A[i][i]
    i=i-1
    #print("ith")

print ("\n Values of a,b,c are: ")
print (x)
print ("\n CPU time : ", time.process_time(),'s')
print ("\n Velocity at t=6s :",x[0]*6*6+x[1]*6+x[2],"m/s")


