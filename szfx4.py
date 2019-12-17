from scipy import integrate
import numpy as np
import math
#def f(x):
    #return math.exp(-x**2)
def f(x):
    if x == 0:
        return 1
    return (math.sin(x))/x
def trapezium(a,b):
    result = ((b-a)/2)*(f(a)+f(b))
    print('%.5f' %result)
    return
def Simpson(a,b):
    result = ((b-a)/6)*(f(a)+4*f((a+b)/2)+f(b))
    print('%.5f' %result)
    return
def Cotes(a,b):
    result = ((b-a)/90)*(7*f(a)+32*f((3*a+b)/4)+12*f((a+b)/2)+32*f((a+3*b)/4)+7*f(b))
    print('%.5f' %result)
    return

def comtrapezium(a,b,e):
    n = 0
    h0 = b-a
    T0 = (f(a)+f(b))*h0/2
    n += 1
    h1 = h0 / 2
    T1 = T0 / 2
    m = int(math.pow(2, n - 1))
    for k in range(1, m + 1):
        T1 += h1 * f(a + (2 * k - 1) * h1)
    while abs(T1-T0) >= e:
        T0 = T1
        h0 = h1
        n += 1
        T1 = T0/2
        h1 = h0/2
        m = int(math.pow(2, n - 1))
        for k in range(1,m+1):
            T1 += h1*f(a+(2*k-1)*h1)
    print('%.6f' % T1)
    print('迭代次数：%d' % n)
    return

def comSimpson(a,b,e):
    n = 0
    h0 = (b-a)/2
    S0 = (f(a)+4*f((a+b)/2)+f(b))*h0/3
    n += 1
    h1 = h0 / 2
    m = int(math.pow(2, n - 1))
    S1 = S0 / 2
    for k in range(1, 2*m + 1):
        S1 += h1 * 4 * f(a + (2 * k - 1) * h1) / 3
    for k in range(1, m + 1):
        S1 -= h1 * 2 * f(a + (4 * k - 2) * h1) / 3
    while abs(S1-S0) >= e:
        S0 = S1
        h0 = h1
        n += 1
        h1 = h0/2
        m = int(math.pow(2, n - 1))
        S1 = S0 / 2
        for k in range(1, 2 * m + 1):
            S1 += h1 * 4 * f(a + (2 * k - 1) * h1) / 3
        for k in range(1, m + 1):
            S1 -= h1 * 2 * f(a + (4 * k - 2) * h1) / 3
    print('%.6f' % S1)
    print('迭代次数：%d' % n)
    return

def Romberg(a,b,e,M):
    T = [0]*M
    S = [0]*M
    h = b - a
    T[0] = h*(f(a)+f(b))/2
    m = 1
    k = 0
    while True:
        h = h/2
        S[0] = 0
        for k in range(1,int(math.pow(2,m-1))+1):
            S[0] += h*f(a+(2*k-1)*h)
        S[0] += T[0] / 2
        k += 1
        for j in range(m):
            S[j+1] = S[j] + (S[j]-T[j])/(math.pow(4,j+1)-1)
        if abs(S[m]-T[m-1]) <= e:
            break
        T = S.copy()
        m += 1
        if m == M:
            print('算法失败')
            return
    print('%.6f' %S[m])
    print('迭代次数：%.d' %k)
    return


#print('梯形求积法：')
#trapezium(0,1)
#print('辛普森求积法：')
#Simpson(0,1)
#print('牛顿——科特斯求积法：')
#Cotes(0,1)
#print('********************')
print('复化梯形公式：')
comtrapezium(0,1,1e-6)
print('********************')
print('复化辛普森公式：')
comSimpson(0,1,1e-6)
#print('********************')
#print('龙贝格算法：')
#Romberg(0,1,1e-5,100)