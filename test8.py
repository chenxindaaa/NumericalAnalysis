import numpy as np


class Test8():


#    def __init__(self):


    def Mifa(self, A, acc, X, N):
        k = 0
        X = X/np.max(X)
        while k < N:
            k += 1
            v = np.dot(A, X)
            res = np.dot(v.transpose(), X)
            v = v/np.max(v)
            if np.max(np.dot(A, v) - np.dot(res, v)) <= acc:
                return res, k
                break
            X = v
        return res, k

if __name__ == '__main__':
    A = np.array([[2, 3, 2], [10, 3, 4], [3, 6, 1]])
    s = Test8()
    print(s.Mifa(A, 1e-5, np.array([1, 1, 1]).reshape(3, 1), 100))