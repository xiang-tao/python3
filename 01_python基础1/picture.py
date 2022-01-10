import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

N = 40
E = 40.0
v = 0.4
a = 0.05
K = 5*10**9
delta = 40*10**(-6)
h = 2*a/N
lam = -2*(1-v**2)/(np.pi*E)
u = 1/K
um = u/h


def C(i, j):
    return (i-j+0.5)*h*(np.log(abs(i-j+0.5)*h)-1)\
           -(i-j-0.5)*h*(np.log(abs(i-j-0.5)*h)-1)


x = np.zeros(N+1)
A = np.zeros((N+1, N+1))
F = np.zeros((N+1, N+1))
b = np.zeros(N+1)
for i in range(N+1):
    b[i] = delta
    F[i][i] = um
    for j in range(N+1):
        A[i][j] = C(x[i], x[j])
# A = A*lam+F
reslut = linalg.solve(A, b)
# reslut[0] = reslut[0]*2
# reslut[N] = reslut[N]*2

plt.plot(x, reslut)
plt.show()