import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches


class Test3:
    def __init__(self):
        fig1 = plt.figure('n = 5')
        ax1 = fig1.add_subplot(111)
        fig2 = plt.figure('n = 10')
        ax2 = fig2.add_subplot(111)
        fig3 = plt.figure('n = 20')
        ax3 = fig3.add_subplot(111)
        self.fig = [ax1, ax2, ax3]

    #  y = 1/(1 + x^2)
    def example(self, x):
        return 1/(1 + x ** 2)

    def example1(self, x):
        return (-2*x) / (1+x**2)**2

    #  拉格朗日插值
    def Lagrange(self, x, y, target):
        res = 0
        for i in range(len(y)):
            l = 1
            for j in x[:i] + x[i+1:]:
                l *= (target - j)/(x[i] - j)
            res += l * y[i]
        return res

    #  画图
    def paint(self, x, y):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        x1 = np.linspace(-5, 5, 100)
        y1 = self.example(x1)
        ax.plot(x1, y1)
        ax.plot(x, y)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

    def Polynomial(self, n):
        old_x = np.linspace(-5, 5, n)
        old_y = self.example(old_x)
        new_x = np.linspace(-5, 5, 100)
        new_y = np.array([])
        for i in new_x:
            new_y = np.append(new_y, self.Lagrange(list(old_x), list(old_y), i))
        self.paint(new_x, new_y)

    #  分段低次多项式插值
    #def piecewise_interpolation(self):

    #  最小二乘法
    def leastSquare(self, x, y):
        if len(x) == 2:
            # 此时x为自然序列
            sx = 0.5 * (x[1] - x[0] + 1) * (x[1] + x[0])
            ex = sx / (x[1] - x[0] + 1)
            sx2 = ((x[1] * (x[1] + 1) * (2 * x[1] + 1))
                   - (x[0] * (x[0] - 1) * (2 * x[0] - 1))) / 6
            x = np.array(range(x[0], x[1] + 1))
        else:
            sx = sum(x)
            ex = sx / len(x)
            sx2 = sum(x ** 2)

        sxy = sum(x * y)
        ey = np.mean(y)

        a = (sxy - ey * sx) / (sx2 - ex * sx)
        b = (ey * sx2 - sxy * ex) / (sx2 - ex * sx)
        return a, b

    # 传入参数格式为np.array,n为阶数
    def leastSquareMulti(self, x, y, n):
        X = [np.sum(x ** i) for i in range(2 * n + 1)]
        Y = np.array([[np.sum(y * x ** i)] for i in range(n + 1)])
        S = np.array([X[i:i + n + 1] for i in range(n + 1)])
        return np.linalg.solve(S, Y)  #

    def triple_corner_method(self, bound, h, target):
        x = bound[0]
        X = x + h
        y = self.example(x)
        Y = self.example(X)
        m = self.example1(x)
        M = self.example1(X)
        res = ((target-X)**2)*(h + 2*(target-x))*y / (h**3) + ((target-x)**2)*(h + 2*(X-target))*Y / (h**3) + \
        ((target-X)**2)*(target-x)*m / (h**2) + ((target-x)**2)*(target-X)*M / (h**2)
        x += h
        X += h
        return res

    def test_1(self):
        x = [11, 12, 13]
        y  = [0.190809, 0.207912, 0.224951]
        print(self.Lagrange(x, y, 11.5))

    def test_2(self):
        for n in [5, 10, 20]:
            self.Polynomial(n)

    def test_3(self):
        n = {5: 'red', 10: 'blue', 20: 'yellow'}
        for i in n:
            count = 0
            h = 10 / i
            for j in range(i):
                target = np.linspace(-5 + count * h, -5 + (count * h) + h, 100)
                res = self.triple_corner_method([-5 + count * h, -5 + (count * h) + h], h, target)
                count += 1
                plt.plot(target, res, color=n[i])
        red_patch = mpatches.Patch(color="red")
        blue_patch = mpatches.Patch(color='blue')
        yellow_patch = mpatches.Patch(color='yellow')
        plt.legend(handles=[red_patch, blue_patch, yellow_patch], labels=["n = 5", 'n = 10', 'n = 20'])
        plt.show()

    def test_4(self):
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([4, 4.5, 6, 8, 8.5])
        a, b = self.leastSquare(x, y)
        plt.scatter(x, y)
        plt.plot(x, a * x + b)
        plt.show()
        E = 0
        for i in range(y.shape[0]):
            E = E + (y[i] - (a*x[i]+b))**2
        print(E)

    def test_5(self):
        x = np.array([2, 3, 4, 7, 8, 10, 11, 14, 16, 18, 19])
        y = np.array([106.42, 108.2, 109.5, 110, 109.94, 110.49, 110.59, 110.6, 110.76, 111, 111.2])
        c = self.leastSquareMulti(x, y, 2)
        plt.scatter(x, y)
        X = np.linspace(0, 20, 200)
        res = c[0] + c[1] * X + c[2] * X**2
        plt.plot(X, res)
        plt.show()
        E = 0
        for i in range(y.shape[0]):
            E = E + (y[i] - (c[0]+c[1]*x[i]+c[2]*x[i]**2)) ** 2
        print(E)

    def test_6(self):
        x = np.array([2, 3, 4, 7, 8, 10, 11, 14, 16, 18, 19])
        y = np.array([106.42, 108.2, 109.5, 110, 109.94, 110.49, 110.59, 110.6, 110.76, 111, 111.2])
        X = 1/x
        Y = np.log(y)
        a, b = self.leastSquare(X, Y)
        plt.scatter(x, y)
        plt.plot(x, np.exp(a * X + b))
        plt.show()
        E = 0
        for i in range(y.shape[0]):
            E = E + (Y[i] - (a*X[i]+b))**2
        print(E)

    def test_7(self):
        n = [5, 10, 20]
        x = np.linspace(-5, 5, 1000)
        y = self.example(x)
        for i in n:
            count = 0
            h = 10/i
            for j in range(i):
                target = np.linspace(-5+count*h, -5+(count*h)+h, 100)
                res = self.triple_corner_method([-5+count*h, -5+(count*h)+h], h, target)
                count += 1
                self.fig[n.index(i)].plot(target, res)
            self.fig[n.index(i)].plot(x, y, '-')
        plt.show()






'''
t.Polynomial(15)
t.leastSquare([1, 2, 3, 4, 5], [4, 4.5, 6, 8, 8.5])
'''
t = Test3()
#t.test_1()
#t.test_2()
#t.test_3()
#t.test_4()
#t.test_5()
#t.test_6()
t.test_7()