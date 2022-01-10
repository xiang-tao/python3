"""
张朝辉cmp流动分析图二代码
"""

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

N = 50
E = 40.0 * 10 ** 6
v = 0.4
a = 0.05
q = 10*10**3
K = 6 * 10 ** 9
# K = 10**11
delta = 60 * 10 ** (-6)
h = 2 * a / N
lam = -2 * (1 - v ** 2) / (np.pi * E)
u = 1 / K
um = u / h


def k(t, s):
    return np.log(abs((t - s) / (-0.025 - s)))


def kk(t, s):
    return np.log(abs((t - s) / (-s)))


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
k1 = np.zeros((N, N))

for i in range(N):
    F[i][i] = r[i]
    b[i] = delta + (1-2*v)*(1+v)/E*q*(x[i]-a)
    for j in range(N):
        if i != j:
            k0[i][j] = k(x[i], x[j])
            k1[i][j] = kk(x[i], x[j])
        if i == j:
            k0[i][j] = (k(x[i] - 0.1 * h, x[j]) + k(x[i] + 0.1 * h, x[j])) / 2
            k1[i][j] = (kk(x[i] - 0.1 * h, x[j]) + kk(x[i] + 0.1 * h, x[j])) / 2

A = h * lam * k0 + um * h * F
A1 = h * lam * k1 + um * h * F

reslut1 = linalg.solve(A, b)
ua = reslut1/K*10**6
up = b*10**6 - ua
plt.plot(x * 10, ua)
plt.plot(x * 10, up)
plt.show()
