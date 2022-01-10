import numpy as np
import scipy.sparse as sp
from scipy import integrate
import scipy.sparse.linalg
import matplotlib.pyplot as plt
from matplotlib import cm
# 组装信息矩阵
a = 0
b = 1
N = 4
h = (b-a)/N
Pb = np.zeros(N+1)
Tb = np.zeros((2, N), dtype=np.int64)
for i in range(N+1):
    Pb[i] = i*h
for i in range(N):
    Tb[0][i] = i
    Tb[1][i] = i+1
P, T = Pb, Tb

print(Tb.dtype)
# 给定基函数极其导数
def d_phi_trail():
    return -1/h


def d_psi_test():
    return 1/h


def c(x):
    return np.e**x


def f(x):
    return d_phi_trail()*d_psi_test()*c(x)


def int(f, x1, x2):
    return (f(x1)+f(x2))*(x2-x1)/2


# 组装矩阵A
A = np.zeros((N+1, N+1))
for i in range(N):
    for j in range(2):
        for k in range(2):
            val = int(f, P[i], P[i+1])
            xi = Tb[k][i]
            yj = Tb[j][i]
            A[xi][yj] += val
print(val)
print(A)










