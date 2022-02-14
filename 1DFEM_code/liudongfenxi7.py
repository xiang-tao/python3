"""
张朝辉cmp流动分析图4代码
"""

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')
import fem

N = 50
E = 40.0 * 10 ** 6
v = 0.4
a = 0.05
K = 2.5 * 10 ** 9
delta = 40 * 10 ** (-6)
q0 = 10*10**3
h = 2 * a / N
lam = -2 * (1 - v ** 2) / (np.pi * E)
u = 1 / K
um = u / h


def k(t, s):
    return np.log(abs((t - s) / (-a - s)))


x = np.zeros(N)
w = np.zeros(N)
r = np.zeros(N)
F = np.zeros((N, N))
b = np.zeros(N)
for i in range(N):
    x[i] = (i + 0.5) * h
    w[i] = 1

r = 1 / w
k0 = np.zeros((N, N))

for i in range(N):
    F[i][i] = r[i]
    b[i] = delta-(1-2*v)*(1+v)/E*q0*(x[i]-a)
    for j in range(N):
        if i != j:
            k0[i][j] = k(x[i], x[j])
        if i == j:
            k0[i][j] = (k(x[i] - 0.5 * h, x[j]) + k(x[i] + 0.5 * h, x[j])) / 2

A = h * lam * k0 + um * h * F
reslut1 = linalg.solve(A, b)
# plt.plot(x * 10, reslut1 / 1000, label='$x_0 = -a$')
# plt.legend()
# plt.show()
ua = reslut1/K*10**6
up = b*10**6 - ua
h0 = 65*10**(-6)
hh = np.zeros(N)
dh = np.zeros(N)
for i in range(N):
    hh[i] = h0*10**6
hh = hh-ua
hh = hh * 10**(-6)
print(hh.shape)
dh[0] = (hh[1]-hh[0])/h
dh[N-1] = (hh[N-1]-hh[N-2])/h
for i in range(1, N-1):
    dh[i] = (hh[i+1]-hh[i-1])/(2*h)

U = 0.4
eta = 0.0025

mesh = fem.IntervalMesh(h, 2*a-h, N-1)
V = fem.FunctionSpace(mesh, degree=1)
u = fem.TrialFunction(V)
v = fem.TestFunction(V)
nabla_u = fem.nabla_u(u)
nabla_v = fem.nabla_v(v)
AA = fem.assemble_matrix(V, nabla_u, nabla_v, hh**3)

ff = -6*eta*U*dh
bb = fem.assemble_vector(V, ff, v)
AA, bb = fem.Dirichlet(0.0, 0.0, AA, bb)
u = fem.Function(V)
fem.solve(AA, u,  bb)
u.generate_value(u.Vector()/80)
yy = u.Vector()
xx = V.generate_pb()
plt.plot(xx, yy)
plt.xlabel("Wafer position(m)", fontsize='15')
plt.ylabel("flow pressure$~p_f$ KPa", fontsize='15')
plt.show()
