import numpy as np


class IntervalMesh(object):
    def __init__(self, left, right, N):
        self.left = left
        self.right = right
        self.N = N

    def generate_p(self):
        P = np.zeros((self.N + 1, 1))
        h = (self.right - self.left) / self.N
        for i in range(self.N + 1):
            P[i][0] = self.left+i * h
        return P

    def generate_t(self):
        T = np.zeros((self.N, 2), dtype=np.int64)
        for i in range(self.N):
            T[i][0] = i
            T[i][1] = i + 1
        return T




