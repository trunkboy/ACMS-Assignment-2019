import numpy

#For matrix system AX=B, solve using gaussion elimination and produce the output along with computational time

A = [[2,1,1,3,2],
[1,2,2,1,1],
[1,2,9,1,5],
[3,1,1,7,1],
[2,1,5,1,8]]

B = [-2,4,3,-5,1]

n=5

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
        print (A)
        print (B)
        print("\n")

        j=j+1

    i=i+1

print (A)
print (B)

#Back Substitution from nth row

x= [0]*n

i = n-1
x[i] = B[i]/A[i][i]
i=i-1
while i>=0:
    sum = 0
    k=i+1
    while k<n:
        print("kth")
        sum = sum + A[i][k]*x[k]
        k=k+1
    x[i]=(B[i]-sum)/A[i][i]
    i=i-1
    print("ith")

print ("answer")
print (x)
