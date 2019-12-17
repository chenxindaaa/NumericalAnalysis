import numpy as np


class Test7:

    def jacobi(self, mat, source, acc):
        for i in range(mat.shape[0]):
            maxi = -mat[i][i]
            mat[i] = mat[i] / maxi
            mat[i][i] = source[i]
            mat[i][-1] *= -1
        k = 0
        while True:
            res = np.dot(mat[:, :-1], source.reshape(source.shape[0], 1)) + mat[:, -1:]
            if np.max(np.abs(res - source)) <= acc:
                break
            source = res
            k += 1
        return res

    def Gauss(self, mat, source, acc):
        A = mat[:, :-1]
        b = mat[:, -1:]
        D = np.diag(np.diag(A))
        U = -1 * np.triu(A, 1)
        L = -1 * (A - np.triu(A, 0))
        temp = np.linalg.inv(D - L)
        k = 0
        while True:
            res = np.dot(np.dot(temp, U), source.reshape(source.shape[0], 1))\
                  + np.dot(temp, b).reshape(source.shape[0], 1)
            if np.max(np.abs(res - source)) <= acc:
                break
            source = res
            k += 1
        return res

    def SOR(self, mat, source, acc, w):
        A = mat[:, :-1]
        b = mat[:, -1:]
        D = np.diag(np.diag(A))
        U = -1 * np.triu(A, 1)
        L = -1 * (A - np.triu(A, 0))
        temp = np.linalg.inv(D - w * L)
        k = 0
        while True:
            res = np.dot(np.dot(temp, (1 - w) * D + w * U), source.reshape(source.shape[0], 1))\
                  + np.dot(w * temp, b).reshape(source.shape[0], 1)
            if np.max(np.abs(res - source)) <= acc:
                break
            source = res
            k += 1
        return k

if __name__ == "__main__":
    mat = np.array([[-4, 1, 1, 1, 1], [1, -4, 1, 1, 1], [1, 1, -4, 1, 1], [1, 1, 1, -4, 1]], dtype=np.float64)
    acc = 1e-5
    source = np.array([0, 0, 0, 0], dtype=np.float64)
    s = Test7()
    for w in np.linspace(1, 10, 100):
        w = 1.2727272727272727
    #print(s.jacobi(mat, source, acc))
    #print(s.Gauss(mat, source, acc))
    #print(s.SOR(mat, source, acc, 1.3))
