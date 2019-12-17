def f(x):
    return x**3+x**2-3*x-3
def dividetwo(a,b,e,g):
    n = 0
    if not f(a) or not f(b):
        x = f(b) if f(a) else f(b)
        print('近似根x=%f\n二分次数n=%d' % (x, n))
        return
    elif f(a)*f(b) < 0:
        while abs(b-a) > e:
            x = (a+b)/2
            if abs(f(x)) < g:
                print ('近似根x=%f\n二分次数n=%d' % (x,n))
                return
            elif f(a)*f(x) < 0:
                b = x
            elif f(b)*f(x) < 0:
                a = x
            n += 1
        x = (a+b)/2
        print('近似根x=%f\n二分次数n=%d' % (x, n))
        return
    return

dividetwo(1.7,1.8,1e-6,1e-9)


