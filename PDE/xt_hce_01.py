import xt_pde_matrix as matrix
import xt_pde_tools as tools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib
"""
向后差分格式求解一维抛物方程代码
"""
e = np.e
pi = np.pi
n = int(input("请输入x轴方向的想要划分的内节点数: "))
times = 3000
tol = 10 ** (-4)
h = 1.0 / (n + 1)
t = int(input("请输入时间轴方向的想要划分的内部层数: "))
tau = 0.1 / (t + 1)
A = np.zeros((n, n))
data = np.zeros(3 * n - 2)
# 考虑到数组不能够用于数组内部，即x[indices[i]]，尽管indices[i]是一个数字，看起来似乎可行，
# 但是当indices[i]是一个数组类型时，将是非法的，indices[i]是列表类型时认为合法
indices = []
indptr = [0]  # 有一个0是因为和算法关系，先初始第一个数
u0 = np.zeros(n)
u1 = np.zeros(n)
# err = np.zeros(times)
B = np.zeros((t + 2, n))
r = tau / (h ** 2)

hx = 0
for i in range(n):
    hx = (i + 1) * h
    u0[i] = e ** (-100 * (hx - 0.25) ** 2) + 0.1 * np.sin(20 * pi * hx)

matrix.back_diff(A, n, r)
matrix.csrmatrix(A, data, indices, indptr)
for i in range(n):
    B[0][i] = u0[i]

for i in range(t + 1):
    tools.csrjacobi(data, indices, indptr, u1, u0, times, tol)
    for k in range(n):
        u0[k] = u1[k]
        B[i + 1][k] = u0[k]

# 画数值解图
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# Make data.
X = np.arange(h, 1, h)
Y = np.arange(0, 0.1 + tau, tau)
X, Y = np.meshgrid(X, Y)
R = X ** 0 + Y ** 0 - 2
Z = R + B

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.rainbow)
fig.colorbar(surf, shrink=0.4, aspect=6)
ax.set_zlim(-1, 1)
ax.contourf(X, Y, Z, zdir='z', offset=-1, cmap=cm.rainbow)
zhfont1 = matplotlib.font_manager.FontProperties\
    (fname="/home/xt/github/python3/python-chinese/SourceHanSansSC-Bold.otf")
# 加标题和轴标签，使用：fontsize=20(any number)可以改变标签字体大小
ax.set_xlabel('X 空间', fontproperties=zhfont1)
ax.set_ylabel('Y 时间', fontproperties=zhfont1)
ax.set_zlabel('Z 温度', fontproperties=zhfont1)
ax.set_title("Title h=1/60 tau=1/600", fontproperties=zhfont1, fontsize=20)
plt.show()
