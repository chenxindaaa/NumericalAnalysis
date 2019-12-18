import numpy as np
def order_Guass_Elimination(A,b,n):
    X = [0] * n
    for i in range(0,n):
        for k in range(i+1,n):
            m = A[k,i]/A[i,i]
            b[0,k] -= m*b[0,i]
            for j in range(i,n):
                A[k,j] -= m*A[i,j]
    p = n-1
    X[p] = b[0,n-1] / A[n-1, n-1]
    p -= 1
    while p >= 0:
        q = p+1
        while q <= n-1:
            X[p] -= A[p,q]*X[q]
            q += 1
        X[p] += b[0,p]
        X[p] /= A[p,p]
        p -= 1
    print('X =', end='')
    print(X)
    return

def List_Guass_Elimination(A,b,n):
    X = [0]*n
    for k in range(n):
        a = max(A[k:n, k])
        if a == 0:
            print('系数矩阵为奇异矩阵')
            return
        for i in range(k,n):
            if A[i,k] == a:
                break
        A[[k, i],k : n] = A[[i, k],k : n]
        b[0,k],b[0,i]= b[0,i],b[0,k]
        for j in range(k+1,n):
            m = A[j, k] / A[k, k]
            b[0,j] -= m * b[0,k]
            for l in range(k,n):
                A[j,l] -= m*A[k,l]
    if A[n-1,n-1] == 0:
        print('系数矩阵为奇异矩阵')
        return
    p = n-1
    X[p] = b[0,n-1] / A[n-1, n-1]
    p -= 1
    while p >= 0:
        q = p+1
        while q <= n-1:
            X[p] -= A[p,q]*X[q]
            q += 1
        X[p] += b[0,p]
        X[p] /= A[p,p]
        p -= 1
    print('X =',end = '')
    print(X)
    return

#A = np.array([[0.101,2.304,3.555],[-1.347,3.712,4.623],[-2.835,1.072,5.643]])
#b = np.array([[1.183,2.137,3.035]])
#A = np.array([[0.3*(1e-15),59.14,3,1],[5.291,-6.13,-1,2],[11.2,9,5,2],[1,2,1,1]])
#b = np.array([[59.17,46.78,1,2]])
A = np.array([[3,2,4],[-1,5,3],[2,6,1]])
b = np.identity(3)
n = 3

print('顺序高斯消去法：')
order_Guass_Elimination(A,b,n)
print('\n')

A = np.array([[0.3*(1e-15),59.14,3,1],[5.291,-6.13,-1,2],[11.2,9,5,2],[1,2,1,1]])
b = np.array([[59.17,46.78,1,2]])
n = 4

print('列主元高斯消去法：')
List_Guass_Elimination(A,b,n)
print('\n')

def LU(A,n,b):
    Ps = np.identity(n)
    Ls = np.identity(n)
    for i in range(n-1):
        P = np.identity(n)
        L = np.identity(n)
        a = max(abs(A[i:n,i]))
        if a == abs(min(A[i:n,i])):
            a *= -1
        for j in range(i,n):
            if A[j,i] == a:
                break
        P[[i,j],:]=P[[j,i],:]
        A = np.dot(P,A)
        for k in range(i+1,n):
            L[k,i] = -A[k,i]/A[i,i]
        A = np.dot(L,A)
        Ps = np.dot(P,Ps)
        Ls = np.dot(Ls,np.dot(P,np.linalg.inv(L)))
    L = np.dot(Ps,Ls)
    U = A
    Y = np.dot(np.linalg.inv(L),np.dot(Ps,b))
    X = np.dot(np.linalg.inv(U),Y)
    result = np.linalg.det(U)*np.linalg.det(Ps)
    print('A的行列式为%d' %result)
    print('A的逆矩阵为：')
    print(X)
    print('L =', end='')
    print(np.float16(L))
    print('U =', end='')
    print(np.float16(U))
    print('线性方程组的解为：')
    print('X=',end='')
    print(X)
    return

print('LU分解：')
LU(A,n,b)

