import xt_pdenew_matrix as matrix
import xt_pde_tools as tools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib

"""
向后差分格式求解数学建模的维抛物方程代码
"""
tem = 65
times = 5000
tol = 10 ** (-6)
t = 3600
tau = 1.0
d = [0.6, 22.9, 3.6, 5.5]  # 每层材料的厚度
n = [6, 8, 6, 6]  # 每层划分的内节点数
h = [0.1 * 10 ** (-3), d[1]/n[1] * 10 ** (-3), 0.6 * 10 ** (-3), d[3]/n[3] * 10 ** (-3)]
# h = 0.2 * 10**(-3)  # 注意此处的单位应该换乘 m
k = [0.082, 0.37, 0.045, 0.028]  # 四层材料对应的热传导率
c = [1377, 2100, 1726, 1005]  # 四层材料了对应的比热容
p = [300, 862, 74.2, 1.18]  # 四层材料对应的密度
r = []
for i in range(4):
    a = k[i] / (c[i] * p[i])
    a = a * (tau / (h[i] ** 2))
    r.append(a)

A1 = matrix.tri_diagonal(n[0], 1 + 2 * r[0], -r[0], -r[0])
A2 = matrix.tri_diagonal(n[1], 1 + 2 * r[1], -r[1], -r[1])
A3 = matrix.tri_diagonal(n[2], 1 + 2 * r[2], -r[2], -r[2])
A4 = matrix.tri_diagonal(n[3], 1 + 2 * r[3], -r[3], -r[3])
A = np.zeros((n[0] + n[1] + n[2] + n[3], n[0] + n[1] + n[2] + n[3]))

A[0:n[0], 0:n[0]] = A1
A[n[0]:n[0] + n[1], n[0]:n[0] + n[1]] = A2
A[n[0] + n[1]:n[0] + n[1] + n[2], n[0] + n[1]:n[0] + n[1] + n[2]] = A3
A[n[0] + n[1] + n[2]:n[0] + n[1] + n[2] + n[3], n[0] + n[1] + n[2]:n[0] + n[1] + n[2] + n[3]] = A4
A[n[0] - 1][n[0]] = -r[0]
A[n[0]][n[0] - 1] = -r[1]
A[n[0] + n[1] - 1][n[0] + n[1]] = -r[1]
A[n[0] + n[1]][n[0] + n[1] - 1] = -r[2]
A[n[0] + n[1] + n[2] - 1][n[0] + n[1] + n[2]] = -r[2]
A[n[0] + n[1] + n[2]][n[0] + n[1] + n[2] - 1] = -r[3]

B = np.zeros((t + 1, n[0] + n[1] + n[2] + n[3]))
B[0] = 37

u0 = np.zeros(n[0] + n[1] + n[2] + n[3])
for i in range(n[0] + n[1] + n[2] + n[3]):
    if i == 0:
        u0[i] = 37 + r[0] * tem
    elif i == n[0] + n[1] + n[2] + n[3] - 1:
        u0[i] = 37 + r[3] * 37
    else:
        u0[i] = 37
u1 = np.zeros(n[0] + n[1] + n[2] + n[3])

data = np.zeros(3 * (n[0] + n[1] + n[2] + n[3]) - 2)
indices = []
indptr = [0]
matrix.csrmatrix(A, data, indices, indptr)
for j in range(t):
    tools.csrjacobi(data, indices, indptr, u1, u0, times, tol)
    B[j + 1] = u1
    for i in range(n[0] + n[1] + n[2] + n[3]):
        if i == 0:
            u0[i] = u1[i] + r[0] * tem
        elif i == n[0] + n[1] + n[2] + n[3] - 1:
            u0[i] = u1[i] + r[3] * (0.077 * 37 + 0.923 * u1[i])
        else:
            u0[i] = u1[i]
print(u1[n[0] + n[1] + n[2] + n[3]-1])
# 画最后时间段温度关于时间的二维曲线图形
x = np.linspace(0, t, t+1)
y = B[:, n[0] + n[1] + n[2] + n[3] - 1]
x1 = np.linspace(0, t-300, t+1-300)
y1 = 44*np.ones(t+1-300)
x2 = 3300*np.ones(t+1-300)
y2 = np.linspace(37, 44, t+1-300)
plt.plot(x, y)
plt.plot(x1, y1)
plt.plot(x2, y2)
zhfont1 = matplotlib.font_manager.FontProperties \
    (fname="/home/xt/PycharmProjects/python_chinese/SourceHanSansSC-Bold.otf")
plt.xlabel("t(时间)", fontproperties=zhfont1)
plt.ylabel("T(温度)", fontproperties=zhfont1)
plt.title("皮肤表面温度", fontproperties=zhfont1)
plt.legend()
plt.show()