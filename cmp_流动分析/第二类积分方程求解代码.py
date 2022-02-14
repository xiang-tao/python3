import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

u = 1.0
lam = 1.0
N = 40
h = 1.0/N
um = u/h
x = np.zeros(N+1)
w = np.zeros(N+1)
r = np.zeros(N+1)
F = np.zeros((N+1,N+1))
b = np.zeros(N+1)
for i in range(N+1):
    x[i] = i*h
    w[i] = 1.0

w[0] = w[0]/2
w[N] = w[N]/2
r = 1/w


def k(t, s):
    return 2*t+s


def g(t):
    return t**2+2/3*t+0.25
    # return 3/2*t+1/3


for i in range(N+1):
    F[i][i] = r[i]
    b[i] = g(x[i])
# print(F)
k0 = np.zeros((N+1, N+1))
for i in range(N+1):
    for j in range(N+1):
        k0[i][j] = k(x[i], x[j])


A = h*lam*k0+um*h*F

reslut = linalg.solve(A, b)
reslut[0] = reslut[0]*2
reslut[N] = reslut[N]*2

plt.plot(x, reslut)
plt.show()