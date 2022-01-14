"""
张朝辉cmp流动分析图4代码
"""

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import frog

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
# A1 = h * lam * k0 + um1 * h * F
reslut1 = linalg.solve(A, b)
# reslut2 = linalg.solve(A1, b)
ua = reslut1/K*10**6
# up = b*10**6 - ua
h0 = 80*10**(-6)
hh = np.zeros(N)
dh = np.zeros(N)
for i in range(N):
    hh[i] = h0*10**6
hh = hh-ua
# plt.plot(x * 10, hh)
# plt.show()
hh = hh * 10**(-6)
dh[0] = (hh[1]-hh[0])/h
dh[N-1] = (hh[N-1]-hh[N-2])/h
for i in range(1,N-1):
    dh[i] = (hh[i+1]-hh[i-1])/(2*h)

U = 0.4
eta = 0.0025

mesh = frog.intervalmesh(h, 2*a-h, N-1)
V = frog.functionspace(mesh, "Lagrange", 1)
u = frog.trialfunction(V)
v = frog.testfunction(V)
nabla_u = frog.nabla_u(u)
nabla_v = frog.nabla_v(v)
a = frog.dot(nabla_u, nabla_v)
AA = frog.assemble_A(a, hh**3)

ff = -6*eta*U*dh
bb = frog.assemble_b(ff, v)
AA, bb = frog.dirichlet1d(0.0, 0.0, AA, bb)
y = frog.solve(AA, bb)
frog.draw(y, V)
plt.show()
