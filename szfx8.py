import numpy as np
def power(A,x,e,N):
    k = 0
    x = x/np.linalg.norm(x)
    while k < N:
        v = np.dot(A,x)
        k += 1
        r = np.dot(v.T,x)
        v = v/np.linalg.norm(v)
        if np.linalg.norm(np.dot(A, v)-r*v) < e:
            break
        x = v
    print('A的按模最大特征值为%.f' %r)
    print('A的按模最大特征值对应的特征向量为:')
    print(np.float16(x))
    print('迭代次数为%d' % k)
    return r,k

def power2(A,x,e,N):
    k = 0
    r = 0
    while k < N:
        x = np.dot(A,x)
        k += 1
        xp = max(abs(x))
        if xp == abs(min(x)):
            r = -xp
            x = -(x / xp)
        else:
            r = xp
            x = x / xp
        if max(abs(np.dot(A, x)-r*x)) < e:
            break
    print(r)
    print('A的按模最大特征值为%.f' %r)
    print('A的按模最大特征值对应的特征向量为:')
    print(np.float16(x))
    print('迭代次数为%d' % k)
    return r,k,x

def repower(A,x,e,N):
    k = 0
    x = x/np.linalg.norm(x)
    r0 = 0
    A1 = np.linalg.inv(A)
    while k < N:
        k += 1
        v = np.dot(A1,x)
        r = 1/np.dot(v.T,x)
        v = v/np.linalg.norm(v)
        if abs(r-r0) <= e:
            break
        x = v
        r0 = r
    print('A的按模最小特征值为%.f' % r)
    print('A的按模最小特征值对应的特征向量为:')
    print(np.float16(x))
    print('迭代次数为%d' % k)
    return r,k

def repower2(A,x,e,N):
    k = 0
    xp = max(abs(x))
    r0 = xp
    A1 = np.linalg.inv(A)
    while k < N:
        k += 1
        x = np.dot(A1,x)
        xp = max(abs(x))
        if xp == abs(min(x)):
            r = -xp
            x = -(x / xp)
        else:
            r = xp
            x = x / xp
        r = 1 / r
        if abs(r-r0) <= e:
            break
        r0 = r
    print('A的按模最小特征值为%.f' % r)
    print('A的按模最小特征值对应的特征向量为:')
    print(np.uint8(x))
    print('迭代次数为%d' % k)

    return r,k,x

def origin_displacement(A,x,e,N,p,k):
    r0 = max(abs(x))
    n = len(x)
    E = np.identity(n)
    A1 = np.linalg.inv(A - p * E)
    while k < N:
        k += 1
        x = np.dot(A1,x)
        xp = max(abs(x))
        if xp == abs(min(x)):
            r = -xp
            x = -(x / xp)
        else:
            r = xp
            x = x / xp
        r = 1 / r
        if abs(r-r0) <= e:
            r += p
            break
        r0 = r
    print('A的按模特征值为%.f' % r)
    print('A的按模特征值对应的特征向量为:')
    print(np.float16((x)))
    print('迭代次数为%d' % k)
    return

def max_Dynamic_origin_displacement(A,x,e1,N1,e2,N2):
    r, k, x = power2(A, x, e1, N1)
    p = r
    origin_displacement(A, x, e2, N2, p,k)
    return

def min_Dynamic_origin_displacement(A,x,e1,N1,e2,N2):
    r, k, x = repower2(A, x, e1, N1)
    p = r
    origin_displacement(A, x, e2, N2, p, k)
    return

A = np.array([[4,-1,1],[-1,3,-2],[1,-2,3]])
x = np.array([[1],[1],[1]])
x = np.array([[2.001],[1.999],[0]])
x = np.array([[0],[1],[1]])
p = 5.9
k = 0


#power(A,x,1e-6,100)
#print('******************')
#power2(A,x,1e-6,100)
#print('******************')
#k = 0
#repower(A,x,1e-6,100)
#print('******************')
#repower2(A,x,1e-6,100)
#print('******************')
#max_Dynamic_origin_displacement(A,x,0.01,10,1e-6,100)
#print('******************')
#min_Dynamic_origin_displacement(A,x,0.1,10,1e-5,100)
#print('******************')
origin_displacement(A,x,1e-6,100,p,k)