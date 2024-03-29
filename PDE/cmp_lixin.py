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

# 请注意r/s单位不是国际单位，应换算成弧度每秒，而一转的弧度是2*pi，故需要乘2*pi
ww = 50.0/60*2*pi  # ----------------------------------------晶片的角速度
wp = 100/60*2*pi  # ----------------------------------------抛光垫的角速度
angle1 = 0.02 * pi / 180  # ---------------------------转角
angle2 = 0.018 * pi / 180  # ---------------------------倾角
d = 0.15  # ---------------------------晶片和抛光垫的旋转中心距
r0 = 5.0*1e-2
p0 = 101000.0
hpiv = 8.0*1e-5  # ----------------------------------------晶片中心高度
viscosity = 0.00214  # ---------------------------抛光液粘度
rho = 1800.0

xx = r0 / hpiv
aa = 6 * viscosity * wp / p0 * xx ** 2
dd = d / r0
ee = ww / wp

nr = 64
n_theta = 64
r = np.linspace(0, 1, nr)
hr = 1 / (nr - 1)
theta = np.linspace(0, 2 * pi, n_theta)
h_theta = 2 * pi / (n_theta - 1)
r_in = r[0:-1]  # 去掉尾
theta_in = theta[0:-1]  # 去掉尾
n_in = (r_in.size - 1) * theta_in.size + 1
D0 = np.zeros(n_in - 1)
D1 = np.zeros(n_in - 1)
Dn = np.zeros(n_in - 1)
D_1 = np.zeros(n_in - 1)
D_n = np.zeros(n_in - 1)
f = np.zeros(n_in)
b = np.zeros(n_in)
A = np.zeros((n_in - 1, n_in - 1))
B = np.zeros((n_in, n_in))
data = np.zeros((theta.size, r.size))
print(data.shape)


def h_function(ri, tj):
    return 1 - xx * r_in[ri] * sin(angle1) * cos(theta_in[tj]) - xx * r_in[ri] \
           * sin(angle2) * sin(theta_in[tj])


def hhh(ri, tj):
    return h_function(ri, tj) ** 3


def dh_r(tj):
    return -(xx * sin(angle1) * cos(theta_in[tj]) + xx * sin(angle2) * sin(theta_in[tj]))


def dh_theta(ri, tj):
    return xx * r_in[ri] * sin(angle1) * sin(theta_in[tj]) - xx * r_in[ri] \
           * sin(angle2) * cos(theta_in[tj])


def dhhh_r(ri, tj):
    return 3 * h_function(ri, tj) ** 2 * dh_r(tj)


def dhhh_theta(ri, tj):
    return 3 * h_function(ri, tj) ** 2 * dh_theta(ri, tj)


def c1(ri, tj):
    return 0.5 * hhh(ri, tj) * hr + r_in[ri] * (hhh(ri, tj) + 0.5 * hr * dhhh_r(ri, tj))


def c2(ri, tj):
    return 2 * r_in[ri] * hhh(ri, tj)


def c3(ri, tj):
    return -0.5 * hhh(ri, tj) * hr + r_in[ri] * (hhh(ri, tj) - 0.5 * hr * dhhh_r(ri, tj))


def c4(ri, tj):
    return hhh(ri, tj) + 0.5 * h_theta * dhhh_theta(ri, tj)


def c5(ri, tj):
    return 2 * hhh(ri, tj)


def c6(ri, tj):
    return hhh(ri, tj) - 0.5 * h_theta * dhhh_theta(ri, tj)


def f1(ri, tj):
    return aa * dd * sin(theta_in[tj]) * r_in[ri] * dh_r(tj)


def f2(ri, tj):
    return aa * (dd * cos(theta_in[tj]) + r_in[ri] + r_in[ri] * ee) * dh_theta(ri, tj)


def lixinf1(ri, tj):
    return 6 * rho * r0 ** 2 / p0 * r_in[ri] ** 2 * h_function(ri, tj) ** 2 \
           * wp ** 2 * dh_r(tj)


def lixinf2(ri, tj):
    return 12 * rho * r0 ** 2 / p0 * r_in[ri] * dd * h_function(ri, tj) ** 2 \
           * cos(theta_in[tj]) * wp ** 2 * dh_r(tj)


def lixinf3(ri, tj):
    return 6 * rho * r0 ** 2 / p0 * dd ** 2 * h_function(ri, tj) ** 2 \
           * cos(theta_in[tj]) ** 2 * wp ** 2 * dh_r(tj)


def lixinf4(ri, tj):
    return -3 * rho * r0 ** 2 / p0 * r_in[ri] ** 2 * h_function(ri, tj) ** 2\
           * ww ** 2 * dh_r(tj)


def lixinf5(ri, tj):
    return -6 * rho * r0 ** 2 / p0 * dd * h_function(ri, tj) ** 2 \
           * sin(theta_in[tj]) * wp ** 2 * dh_theta(ri, tj)


def lixinf6(ri, tj):
    return -3 * rho * r0 ** 2 / p0 / r_in[ri] * dd ** 2 * h_function(ri, tj) ** 2 \
           * sin(2 * theta_in[tj]) * wp ** 2 * dh_theta(ri, tj)


def lixinf7(ri, tj):
    return 4 * rho * r0 ** 2 / p0 * r_in[ri] * hhh(ri, tj) * wp ** 2


def lixinf8(ri, tj):
    return -2 * rho * r0 ** 2 / p0 * r_in[ri] * hhh(ri, tj) * ww ** 2


def lixinf9(ri, tj):
    return 2 * rho * r0 ** 2 / p0 * dd * hhh(ri, tj) * cos(theta_in[tj]) * wp ** 2


def lixinf10(ri, tj):
    return -2 * rho * r0 ** 2 / (p0 * r_in[ri]) * wp ** 2 * hhh(ri, tj) \
           * dd ** 2 * cos(2 * theta_in[tj])


def lixinf(ri, tj):
    return lixinf1(ri, tj) + lixinf2(ri, tj) + lixinf3(ri, tj) + lixinf4(ri, tj)\
           + lixinf5(ri, tj) + lixinf6(ri, tj) + lixinf7(ri, tj) + lixinf8(ri, tj) \
           + lixinf9(ri, tj) + lixinf10(ri, tj)


def suma():
    summ = 0
    for j in range(len(theta_in)):
        summ += hhh(1, j)
    return summ


def sumb():
    summ = 0
    for j in range(len(theta_in)):
        summ += sin(theta_in[j]) * h_function(1, j)
    return -aa * dd * hr * summ


k = 0
for i in range(1, r_in.size):
    for j in range(theta_in.size):
        D0[k] = c2(i, j) / hr ** 2 + c5(i, j) / r_in[i] / h_theta ** 2
        D1[k] = -c4(i, j) / r_in[i] / h_theta ** 2
        D_1[k] = -c6(i, j) / r_in[i] / h_theta ** 2
        Dn[k] = -c1(i, j) / hr ** 2
        D_n[k] = -c3(i, j) / hr ** 2
        b[k + 1] = -f1(i, j) - f2(i, j) - lixinf(i, j)
        k += 1

b[0] = sumb()

for i in range(n_in - len(theta_in), n_in):
    f[i] = Dn[i - 1]

matrix.fill_diag(A, D0, 0)
matrix.fill_diag(A, D1, 1)
matrix.fill_diag(A, D_1, -1)
matrix.fill_diag(A, D_n, -theta_in.size)
matrix.fill_diag(A, Dn, theta_in.size)

for i in range(n_in - 1):
    if i % len(theta_in) == 0:
        if i == 0:
            A[i][i + len(theta_in) - 1] = D_1[i]
        else:
            A[i][i + len(theta_in) - 1] = A[i][i - 1]
            A[i][i - 1] = 0
    if (i + 1) % len(theta_in) == 0:
        if i == n_in - 2:
            A[i][i - len(theta_in) + 1] = D1[i]
        else:
            A[i][i - len(theta_in + 1) + 1] = A[i][i + 1]
            A[i][i + 1] = 0

B[0][0] = suma()
for j in range(len(theta_in)):
    B[0][j + 1] = -hhh(1, j)
    B[j + 1][0] = D_n[j]
B[1:B.shape[0], 1:B.shape[1]] = A

rb = b - f
rx = linalg.solve(B, rb)

k = 1
for j in range(1, r.size - 1):
    for i in range(theta.size - 1):
        data[i][j] = rx[k]
        k += 1
data[:, 0] = rx[0]
data[:, data.shape[1] - 1] = 1.0
data[data.shape[0] - 1, :] = data[0, :]

print(data.max())
print(data.min())

X, Y = np.meshgrid(r, theta)

fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
surf = ax.contourf(Y, X, data, cmap=cm.rainbow)
fig.colorbar(surf)

X1, Y1 = np.meshgrid(np.linspace(0, 1, nr + 1), np.linspace(0, 2 * pi, n_theta + 1))
xx1 = X1 * np.cos(Y1)
yy1 = X1 * np.sin(Y1)
fig1, ax1 = plt.subplots()
surf1 = ax1.pcolor(xx1, yy1, data, cmap=cm.rainbow)
fig1.colorbar(surf1)

xx = X * np.cos(Y)
yy = X * np.sin(Y)
fig2, ax2 = plt.subplots(subplot_kw={"projection": "3d"})
surf2 = ax2.plot_surface(xx, yy, data, cmap=cm.rainbow)

xx = xx.transpose()
xx = xx.reshape(xx.shape[0] * xx.shape[1], )

yy = yy.transpose()
yy = yy.reshape(yy.shape[0] * yy.shape[1], )

zz = data.transpose()
zz = zz.reshape(data.shape[0] * data.shape[1], )

data1 = np.array([xx, yy, zz])

data1 = data1.transpose()

np.savetxt('/home/xt/github/python3/file_and_animation/cmpdata_lixin64.plt', np.c_[data1],
           fmt='%.16f', delimiter='\t')

plt.show()
