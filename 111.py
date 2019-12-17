import matplotlib.pyplot as plt
import numpy as np
import math
from numpy import matlib
import scipy.stats as stats


class a():
    def __init__(self):
        self.fig = plt.figure()
        self.axes = self.fig.add_subplot(111)


    def start_Single_Normal_plot1(self, m, n, mu, sigma):
        self.axes.cla()
        X = sigma * matlib.randn((n, m)) + mu  # 模拟抽样，每一列是一个样本容量为n的样本，共抽样m次
        Xmean = np.mean(X, axis=0)  # 样本均值
        z = (Xmean - mu) / (sigma / math.sqrt(n))
        Z = np.array(z).flatten()
        x = np.linspace(-5, 5, 100)
        N, bins, patches = self.axes.hist(Z, 50, density=1, facecolor="yellow", edgecolor="black", alpha=0.7)
        self.axes.plot(x, stats.norm.pdf(x), "b", label="N(0,1)")
        self.axes.set_xlabel('x')
        self.axes.set_ylabel('pdf')
        self.fig.suptitle("方差$\sigma^2$已知时,样本均值的抽样分布")
        self.axes.legend()
        plt.show()


if __name__ == '__main__':
    a = a()
    a.start_Single_Normal_plot1(10000, 20, 0, 1)