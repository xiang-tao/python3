"""
张朝辉cmp流动分析图4代码
"""

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

N = 50
E = 40.0 * 10 ** 6
v = 0.4
a = 0.05
K = 1 * 10 ** 9
K1 = 1 * 10 ** 10
delta = 60 * 10 ** (-6)
h = 2 * a / N
lam = -2 * (1 - v ** 2) / (np.pi * E)
u = 1 / K
u1 = 1/K1
um = u / h
um1 = u1/h


def k(t, s):
    return np.log(abs((t - s) / (-0.025 - s)))


x = np.zeros(N)
w = np.zeros(N)
r = np.zeros(N)
F = np.zeros((N, N))
b = np.zeros(N)
for i in range(N):
    # x[i] = (i+0.001)*h
    x[i] = (i + 0.5) * h
    w[i] = 1

r = 1 / w
k0 = np.zeros((N, N))

for i in range(N):
    F[i][i] = r[i]
    b[i] = delta
    for j in range(N):
        if i != j:
            k0[i][j] = k(x[i], x[j])
        if i == j:
            k0[i][j] = (k(x[i] - 0.1 * h, x[j]) + k(x[i] + 0.1 * h, x[j])) / 2

A = h * lam * k0 + um * h * F
A1 = h * lam * k0 + um1 * h * F
reslut1 = linalg.solve(A, b)
reslut2 = linalg.solve(A1, b)
ua = reslut2/K1*10**6
up = b*10**6 - ua
print(up)
plt.plot(x * 10, ua)
plt.plot(x * 10, up)
plt.show()
