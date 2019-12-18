import numpy as np


class Test7:

    def jacobi(self, mat, source, acc):
        for i in range(mat.shape[0]):
            maxi = -mat[i][i]
            mat[i] = mat[i] / maxi
            mat[i][i] = 0
            mat[i][-1] *= -1
        k = 0
        while True:
            res = np.dot(mat[:, :-1], source.reshape(source.shape[0], 1)) + mat[:, -1:]
            if np.max(np.abs(res - source)) <= acc:
                break
            source = res
            k += 1
        return res, k

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
        return res, k

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

    def test_1(self):
        for i in [[0, 0, 0], [1, 1, 1],  [2, 2, 2]]:
            acc = 1e-5
            mat = np.array([[5, 2, 1, 8], [2, 8, -3, 21], [1, -3, -6, 1]], dtype=np.float64)
            source = np.array(i, dtype=np.float64)
            print('当初始向量为', i, '时：')
            res, k = s.jacobi(mat, source, acc)
            print(res)
            print('迭代次数为%d' % k)

    def test_2(self):
        for i in [[0, 0, 0], [1, 1, 1],  [2, 2, 2]]:
            acc = 1e-5
            mat = np.array([[5, 2, 1, 8], [2, 8, -3, 21], [1, -3, -6, 1]], dtype=np.float64)
            source = np.array(i, dtype=np.float64)
            print('当初始向量为', i, '时：')
            res, k = s.Gauss(mat, source, acc)
            print(res)
            print('迭代次数为%d' % k)

    def test_3(self):
        mat = np.array([[-4, 1, 1, 1, 1], [1, -4, 1, 1, 1], [1, 1, -4, 1, 1], [1, 1, 1, -4, 1]], dtype=np.float64)
        acc = 1e-5
        source = np.array([0, 0, 0, 0], dtype=np.float64)
        res, k = s.jacobi(mat, source, acc)
        print('雅可比迭代法：')
        print(res)
        print('迭代次数： %d ' % k)
        mat = np.array([[-4, 1, 1, 1, 1], [1, -4, 1, 1, 1], [1, 1, -4, 1, 1], [1, 1, 1, -4, 1]], dtype=np.float64)
        acc = 1e-5
        source = np.array([0, 0, 0, 0], dtype=np.float64)
        res, k = s.Gauss(mat, source, acc)
        print('高斯赛德尔迭代法：')
        print(res)
        print('迭代次数： %d ' % k)

    def test_4(self):
        for i in range(1, 20):
            mat = np.array([[-4, 1, 1, 1, 1], [1, -4, 1, 1, 1], [1, 1, -4, 1, 1], [1, 1, 1, -4, 1]], dtype=np.float64)
            acc = 1e-5
            source = np.array([0, 0, 0, 0], dtype=np.float64)
            i /= 10
            print(s.SOR(mat, source, acc, i), '  w = %f' % i)

    def test_5(self):
        for i in range(1, 15):
            i /= 10
            acc = 1e-5
            mat = np.array([[1, 2, -2, 8], [2, 8, -3, 21], [1, -3, -6, 1]], dtype=np.float64)
            source = np.array([0, 0, 0], dtype=np.float64)
            print(s.SOR(mat, source, acc, i), '\tw = %f' % i)

    def test_6(self):
        acc = 1e-5
        mat = np.array([[2, -1, 1, 2], [2, 2, 2, 6], [-1, -1, 2, 0]], dtype=np.float64)
        source = np.array([0, 0, 0], dtype=np.float64)
        res, k = s.Gauss(mat, source, acc)
        print('高斯赛德尔：', res,  'k = %d' % k)






if __name__ == "__main__":
    s = Test7()
    '''
    mat = np.array([[-4, 1, 1, 1, 1], [1, -4, 1, 1, 1], [1, 1, -4, 1, 1], [1, 1, 1, -4, 1]], dtype=np.float64)
    acc = 1e-5
    source = np.array([0, 0, 0, 0], dtype=np.float64)
    for w in np.linspace(1, 10, 100):
        w = 1.2727272727272727
    '''
    #print(s.jacobi(mat, source, acc))
    #print(s.Gauss(mat, source, acc))
    #print(s.SOR(mat, source, acc, 1.3))
    s.test_6()
