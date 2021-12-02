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
beta = np.linspace(-0.03, 0.03, 7)
npwf = np.zeros(7)
npwx = np.zeros(7)
npwy = np.zeros(7)

# 请注意r/s单位不是国际单位，应换算成弧度每秒，而一转的弧度是2*pi，故需要乘2*pi
ww = 50.0 / 60 * 2 * pi  # ----------------------------------------晶片的角速度
wp = 100.0 / 60 * 2 * pi  # ----------------------------------------抛光垫的角速度
angle1 = 0.02 * pi / 180  # ---------------------------转角
# angle2 = 0.018 * pi / 180  # ---------------------------倾角
d = 0.15  # ---------------------------晶片和抛光垫的旋转中心距
r0 = 5.0 * 1e-2
p0 = 101000.0
hpiv = 8.0 * 1e-5  # ----------------------------------------晶片中心高度
viscosity = 0.00214  # ---------------------------抛光液粘度
rho = 1800.0

for num in range(7):
    angle2 = beta[num]*pi/180
    xx = r0 / hpiv
    aa = 6 * viscosity * wp / p0 * xx ** 2
    dd = d / r0
    ee = ww / wp

    nr = 32
    n_theta = 32
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


    def duiliuf1(ri, tj):
        return 6 * rho * r0 ** 2 * wp ** 2 / p0 * r_in[ri] * (r_in[ri] + dd
                                                              * cos(theta_in[tj])) * h_function(ri, tj) ** 2 * dh_r(tj)


    def duiliuf2(ri, tj):
        return -3 * ee ** 2 * rho * r0 ** 2 * wp ** 2 / p0 * r_in[ri] ** 2 \
               * h_function(ri, tj) ** 2 * dh_r(tj)


    def duiliuf3(ri, tj):
        return -6 * dd * rho * r0 ** 2 * wp ** 2 / p0 * h_function(ri, tj) ** 2 \
               * sin(theta_in[tj]) * dh_theta(ri, tj)


    def duiliuf4(ri, tj):
        return 2 * rho * r0 ** 2 * wp ** 2 * hhh(ri, tj) / p0 * (2 * r_in[ri] + dd * cos(theta_in[tj]))


    def duiliuf5(ri, tj):
        return -2 * ee ** 2 * rho * r0 ** 2 * wp ** 2 * r_in[ri] * hhh(ri, tj) / p0


    def duiliuf6(ri, tj):
        return -2 * dd * rho * r0 ** 2 * wp ** 2 / p0 * hhh(ri, tj) * cos(theta_in[tj])


    def duiliuf(ri, tj):
        return duiliuf1(ri, tj) + duiliuf2(ri, tj) + duiliuf3(ri, tj) \
               + duiliuf4(ri, tj) + duiliuf5(ri, tj) + duiliuf6(ri, tj)


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
            b[k + 1] = -f1(i, j) - f2(i, j) - duiliuf(i, j)
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

    wf = 0.0
    wx = 0.0
    wy = 0.0
    for j in range(n_theta):
        for i in range(nr):
            wf += data[j][i] * r[i] * h_theta * hr
            wx += data[j][i] * r[i] ** 2 * sin(theta[j]) * h_theta * hr
            wy += data[j][i] * r[i] ** 2 * cos(theta[j]) * h_theta * hr
    npwf[num] = wf / pi
    npwx[num] = wx / pi
    npwy[num] = -wy / pi

cmpwww = np.array([beta, npwx, npwy, npwf])

np.savetxt('/home/xt/github/python3/cmpdata/beta_www32.txt', np.c_[cmpwww],
           fmt='%.16f', delimiter='\t')
