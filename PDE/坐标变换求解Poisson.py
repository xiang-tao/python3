import xt_pdenew_matrix as matrix
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from matplotlib import cm

"""
利用坐标变换，在极坐标下求解原方程解为u = x^2+y^2-1的拉普拉斯方程，即-delta u = -4
"""

sin = np.sin
cos = np.cos
pi = np.pi

nr = 32
n_theta = 32
r = np.linspace(0, 1, nr)
hr = 1 / (r.size - 1)
theta = np.linspace(0, 2 * pi, n_theta)
h_theta = 2 * pi / (theta.size - 1)
r_in = r[1:-1]  # 去掉首尾
theta_in = theta[1:-1]  # 去掉首尾
n_in = r_in.size * theta_in.size
D0 = np.zeros(n_in)
D1 = np.zeros(n_in)
Dn = np.zeros(n_in)
D_1 = np.zeros(n_in)
D_n = np.zeros(n_in)
f = np.zeros(n_in)
b = np.zeros(n_in)
A = np.zeros((n_in, n_in))

k = 0
for i in range(n_in):
    D0[i] = 2 / hr ** 2 + 2 / r_in[k] ** 2 / h_theta ** 2
    D1[i] = -1 / r_in[k] ** 2 / h_theta ** 2
    D_1[i] = -1 / r_in[k] ** 2 / h_theta ** 2
    Dn[i] = -(r_in[k] + hr / 2) / r_in[k] / hr ** 2
    D_n[i] = -(r_in[k] - hr / 2) / r_in[k] / hr ** 2
    b[i] = -4
    if (i+1) % theta_in.size == 0:
        k += 1

# D1
k = theta_in.size-1
for i in range(r_in.size):
    f[k] += -D1[k]*(r_in[i]**2-1)
    D1[k] = 0
    k += theta_in.size

# D_1
k = 0
for i in range(r_in.size):
    f[k] += (-D_1[k] * (r_in[i] ** 2 - 1))
    D_1[k] = 0
    k += theta_in.size

# Dn
k = (r_in.size-1)*theta_in.size
for i in range(theta_in.size):
    # f[k] += -Dn[k]*0
    Dn[k] = 0
    k += 1

# D_n
k = 0
for i in range(theta_in.size):
    f[k] += -D_n[k]*(-1)
    D_n[k] = 0
    k += 1

matrix.fill_diag(A, D0, 0)
matrix.fill_diag(A, D1, 1)
matrix.fill_diag(A, D_1, -1)
matrix.fill_diag(A, D_n, -theta_in.size)
matrix.fill_diag(A, Dn, theta_in.size)

rb = b+f
rx = linalg.solve(A, rb)

B = np.zeros((theta_in.size + 2, r_in.size + 2))
B[:, 0:1] = -1
k = 0
for i in range(r_in.size):
    for j in range(theta_in.size):
        B[j + 1][i + 1] = rx[k]
        k += 1
B[0:1, 1:r_in.size + 1] = r_in ** 2 - 1
B[theta_in.size + 1:theta_in.size + 2, 1:r_in.size + 1] = r_in ** 2 - 1

X, Y = np.meshgrid(r, theta)

fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
surf = ax.contourf(Y, X, B, cmap=cm.rainbow)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

xx = X * np.cos(Y)
yy = X * np.sin(Y)

xx = xx.transpose()
xx = xx.reshape(xx.shape[0] * xx.shape[1], )

yy = yy.transpose()
yy = yy.reshape(yy.shape[0] * yy.shape[1], )

zz = B.transpose()
zz = zz.reshape(B.shape[0] * B.shape[1], )

data = np.array([xx, yy, zz])

data = data.transpose()


# np.savetxt('/home/xt/github/python3/file_and_animation/matrix.plt', np.c_[data],
#            fmt='%.8f', delimiter='\t')
