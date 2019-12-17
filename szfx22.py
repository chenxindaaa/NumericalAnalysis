def f(x):
    return x**3+x**2-3*x-3
def df(x):
    return 3*(x**2)+2*x-3
def newton(x,e,g,N):
    if not df(x):
        print('算法失败')
        return
    n = 0
    while n <= N:
        x1 = x - f(x)/df(x)
        n += 1
        if abs(x1-x) < e or abs(f(x1)) < g:
            print('牛顿迭代法:\n近似根x=%f\n迭代次数n=%d' % (x1, n))
            return
        x = x1
    print('超过最大迭代次数')
    return
def singlenewton(x0,x1,e,g,N):
    if not df(x0) or not df(x1):
        print('算法失败')
        return
    n = 0
    while n <= N:
        x2 = x1 - f(x1)*(x1-x0)/(f(x1)-f(x0))
        n += 1
        if abs(x2-x1) < e or abs(f(x2)) < g:
            print('单点弦截法:\n近似根x=%f\n迭代次数n=%d' % (x2, n))
            return
        x1 = x2
    print('超过最大迭代次数')
    return
def doublenewton(x0,x1,e,g,N):
    if not df(x0) or not df(x1):
        print('算法失败')
        return
    n = 0
    while n <= N:
        x2 = x1 - f(x1)*(x1-x0)/(f(x1)-f(x0))
        n += 1
        if abs(x2-x1) < e or abs(f(x2)) < g:
            print('两点弦截法:\n近似根x=%f\n迭代次数n=%d' % (x2, n))
            return
        x0 = x1
        x1 = x2
    print('超过最大迭代次数')
    return

newton(1.5,1e-6,1e-9,30)
print('********************')
singlenewton(1.5,1.7,1e-6,1e-9,50)
print('********************')
doublenewton(1.5,1.7,1e-6,1e-9,30)