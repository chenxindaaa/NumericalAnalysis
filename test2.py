import numpy as np
from matplotlib import pyplot as plt


class Test1:
    def __init__(self, fun):
        if fun == 'fun1':
            self.fun = lambda x: x ** 3 + x ** 2 - 3
            self.fun1 = lambda x: 3 * x ** 2 + 2 * x
            self.fig = plt.figure()
            self.axes = self.fig.add_subplot(111)


    def double_division(self, x, acc1, acc2):
        n = 0
        mini = self.fun(x[0])
        maxi = self.fun(x[1])
        if mini * maxi > 0:
            print('算法失效')
        else:
            while abs(maxi - mini) > acc1:
                res = (x[0] + x[1])/2
                f = self.fun(res)
                if abs(f) < acc2:
                    return res, n
                    break
                if f * maxi < 0:
                    x[0] = res
                    mini = f
                else:
                    x[1] = res
                    maxi = f
                n += 1
        return (x[0] + x[1])/2, n


    def Newton(self, x, acc1, acc2, N):
        n = 0
        if self.fun1(x) == 0 or n > N:
            print('算法失效')
        else:
            temp = float("inf")
            while not abs(x - temp) < acc1 or not abs(self.fun(x)) < acc2:
                temp = x
                x = x - self.fun(x)/self.fun1(x)
                n += 1
            return x, n


    def one_point(self, x0, x, acc1, acc2):
        n = 0
        f0 = self.fun(x0)
        if self.fun(x) == 0:
            print('算法失效')
        else:
            temp = float("inf")
            while not abs(x - temp) < acc1 or not abs(self.fun(x)) < acc2:
                temp = x
                x = x - (self.fun(x) / (self.fun(x) - f0)) * (x - x0)
                n += 1
            return x, n

    def two_point(self, x0, x, acc1, acc2):
        n = 0
        if self.fun(x) == 0:
            print('算法失效')
        else:
            while not abs(x - x0) < acc1 or not abs(self.fun(x)) < acc2:
                temp = x
                x = x - (self.fun(x) / (self.fun(x) - self.fun(x0))) * (x - x0)
                x0 = temp
                n += 1
            return x, n


if __name__ == "__main__":
    s = Test1('fun1')
    print(s.double_division([1, 3], 1e-6, 1e-9))
    #print(s.Newton(1.5, 1e-6, 1e-9, 100))
    #print(s.one_point(3, 2, 1e-6, 1e-9))
    #print(s.two_point(3, 2, 1e-6, 1e-9))



