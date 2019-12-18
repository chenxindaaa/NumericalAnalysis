import numpy as np
import matplotlib.pyplot as plt


class Test5:

    def __init__(self, fun):
        if fun == 'fun1':
            self.fun = lambda x, y:y ** 2
        elif fun == 'fun2':
            self.fun = lambda x, y:x / y
        elif fun == 'fun3':
            self.fun = lambda x, y: -50 * y
        elif fun == 'fun4':
            self.fun = lambda x, y: np.sin(x) / x

    def paint(self, x, y):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        x0 = np.linspace(0, 0.4, 50)
        y0 = 1 / (1 - x0)
        print(x)
        print(y)
        print(1 / (1 - x))
        ax.plot(x0, y0, label='origin')
        ax.plot(x, y, label='RungeKutta')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()
        plt.show()

    def ObviousEuler(self,xn, X, yn, h):
        y = [yn]
        while xn <= X :
            yn = yn + h * self.fun(xn, yn)
            y.append(yn)
            xn += h
        return y

    def hidingEuler(self, xn, X, yn, h):
        y = [yn]
        while xn <= X :
            yn = yn / (1 + 50 * h)
            y.append(yn)
            xn += h
        return y


    def improvedEuler(self,xn, X, yn, h):
        y = [yn]
        while xn < X:
            Y = yn + h * self.fun(xn, yn)
            yn = yn + h / 2 * (self.fun(xn, yn) + self.fun(xn + h, Y))
            y.append(yn)
            xn += h
        return y

    def RKF(self,xn, X, yn, h):
        y = [yn]
        while xn < X:
            F1 = h * self.fun(xn, yn)
            F2 = h * self.fun(xn + h / 4, yn + F1 / 4)
            F3 = h * self.fun(xn + 3 / 8 * h, yn + 3 / 32 * F1 + 9/32*F2)
            F4 = h * self.fun(xn + 12/13*h, yn+1932/2197*F1-7200/2197*F2+7296/2197*F3)
            F5 = h * self.fun(xn + h, yn + 439/216*F1-8*F2+3680/513*F3-845/4104*F4)
            F6 = h * self.fun(xn + 0.5*h, yn - 8/27*F1+2*F2-3544/2565*F3+1859/4104*F4-11/40*F5)
            yn = yn + 16/135*F1 + 6656/12825*F3 + 28561/56430*F4 - 9/50*F5 +2/55*F6
            y.append(yn)
            xn += h
        return y

    def RungeKutta(self,xn, X, yn, h):
        y = [yn]
        while xn < X:
            k1 = self.fun(xn, yn)
            k2 = self.fun(xn + 1 / 2 * h, yn + 1 / 2 * h * k1)
            k3 = self.fun(xn + 1 / 2 * h, yn + 1 / 2 * h * k2)
            k4 = self.fun(xn + h, yn + h * k3)
            yn = yn + h/6 * (k1 + 2*k2 + 2*k3 + k4)
            y.append(yn)
            xn += h
        return y

    def test_1(self):
        y = s.RungeKutta(0, 0.4, 1, 0.1)
        #y = s.improvedEuler(0, 0.4, 1, 0.1)
        x = np.linspace(0, 0.4, 5)
        s.paint(x, y)

    def test_2(self):
        y = s.ObviousEuler(0, 0.02, 100, 0.001)
        x = np.linspace(0, 0.02, 21)
        y0 = 100 * np.e ** (-50*x)
        print(x)
        print(y)
        print(y0)
        plt.scatter(x, y0, label='origin')
        plt.plot(x, y, label='ObviousEuler h=0.001')
        plt.legend()
        plt.show()

    def test_3(self):
        for h in [0.01, 0.005, 0.001]:
            y = s.ObviousEuler(0.00000000000000000000000000000001, 1, 1, h)
            print('h = %f 时，近似值为 %.4f' % (h, (y[-1] - y[0])))

    def test_4(self):
        y = s.RKF(0, 0.4, 1, 0.1)
        x = [0.1*i for i in range(5)]
        for j in range(len(x)):
            print('x: %f,y: %f' % (x[j], y[j]))

    def test_5(self):
        y = s.RKF(2, 2.6, 1, 0.1)
        x = [0.1 * i for i in range(20, 27)]
        for j in range(len(x)):
            print('x: %f,y: %f' % (x[j], y[j]))




if __name__ == '__main__':
    s = Test5('fun2')
    s.test_5()

