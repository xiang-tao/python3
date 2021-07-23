import xt_pdenew_matrix as matrix
import xt_pde_tools as tools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib

"""
向后差分格式求解数学建模的维抛物方程代码
"""
times = 5000
tol = 10 ** (-4)
t = 5400
tau = 1.0
d = [0.6, 6, 3.6, 5]  # 每层材料的厚度
n = [6, 10, 6, 10]  # 每层划分的内节点数
h = [0.1 * 10 ** (-3), 0.6 * 10 ** (-3), 0.6 * 10 ** (-3), 0.5 * 10 ** (-3)]
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
        u0[i] = 37 + r[0] * 75
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
            u0[i] = u1[i] + r[0] * 75
        elif i == n[0] + n[1] + n[2] + n[3] - 1:
            u0[i] = u1[i] + r[3] * (0.077 * 37 + 0.923 * u1[i])
        else:
            u0[i] = u1[i]

"""
# 画最后时间段温度关于时间的二维曲线图形
x = np.linspace(0, t, t+1)
y = B[:, n[0] + n[1] + n[2] + n[3] - 1]
plt.plot(x, y)
plt.show()

"""

# 画三维图形
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# Make data.
# X = np.arange(h, 15.4, h)
X = np.zeros(32)
k = 0
for i in range(n[0]):
    X[k] = (i + 1) * h[0] * 10 ** 3
    k += 1
for i in range(n[1]):
    X[k] = 0.6 + (i + 1) * h[1] * 10 ** 3
    k += 1
for i in range(n[2]):
    X[k] = 6.6 + (i + 1) * h[2] * 10 ** 3
    k += 1
for i in range(n[3]):
    X[k] = 10.2 + (i + 1) * h[3] * 10 ** 3
    k += 1
Y = np.arange(0, t + 1, tau)
X, Y = np.meshgrid(X, Y)
R = X ** 0 + Y ** 0 - 2
Z = R + B

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.rainbow)
fig.colorbar(surf, shrink=0.4, aspect=6)

# 使用中文
zhfont1 = matplotlib.font_manager.FontProperties \
    (fname="/home/xt/PycharmProjects/python_chinese/SourceHanSansSC-Bold.otf")

# 加标题和轴标签
ax.set_xlabel('X 轴(空间)', fontproperties=zhfont1)
ax.set_ylabel('Y 轴(时间)', fontproperties=zhfont1)
ax.set_zlabel('Z 温度', fontproperties=zhfont1)
ax.set_title("衣服外面至皮肤表面温度分布数值解", fontproperties=zhfont1, fontsize=20)

plt.show()


