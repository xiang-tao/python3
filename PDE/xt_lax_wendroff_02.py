import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib

"""
当a<0时候的lax-wendroff格式，左边界取u(k)_0=u(k)_1
"""

pi = np.pi
sin = np.sin
d = int(input("请输入空间方向的长度: "))
t = int(input("请输入想要终止的时间: "))
nx = int(input("请输入空间方向的想要划分的节点数: "))
nt = int(input("请输入时间方向的想要划分的节点数: "))
a = int(input("请输入系数a: "))
tau = t / (nt - 1)
h = d / (nx - 1)
r = a * tau / h
u0 = np.zeros(nx - 2)
ut = np.zeros(nt)
B = np.zeros((nt, nx))
A = np.zeros((nx - 2, nx - 2))
for i in range(nx - 2):
    if i == 0:
        A[i][i] = 1 - r * r
        A[i][i + 1] = r * (r - 1) / 2
    if i == nx - 3:
        A[i][i] = 1 - r * r
        A[i][i - 1] = r * (r + 1) / 2
    else:
        A[i][i] = 1 - r * r
        A[i][i - 1] = r * (r + 1) / 2
        A[i][i + 1] = r * (r - 1) / 2
for i in range(nt):
    mt = i * tau
    ut[i] = 1 + sin(4 * pi * mt)
for i in range(nx - 2):
    mx = (i + 1) * h
    u0[i] = 1 + sin(2 * pi * mx)
B[:, nx - 1] = ut
B[0, 1: nx - 1] = u0
u2 = [1]  # 左边界值的第一个初值点为１

for i in range(nt - 1):
    u1 = A @ u0
    u1[0] += r * (r + 1) * u2[i] / 2
    u1[nx - 3] += r * (r - 1) * ut[i] / 2
    k = u1[0]
    u2.append(k)
    u0 = u1
    B[i + 1, 1: nx - 1] = u0

B[:, 0] = u2
# 获取真解数据
C = np.zeros((nt, nx))
for i in range(nt):
    zht = i * tau
    for j in range(nx):
        zhx = j * h
        C[i][j] = 1 + sin(2 * pi * (zhx + 2 * zht))

fig = plt.figure(figsize=plt.figaspect(0.5))
# Make data.
X = np.arange(0, d + h, h)
Y = np.arange(0, t + tau, tau)
X, Y = np.meshgrid(X, Y)
R = X ** 0 + Y ** 0 - 2
Z = R + B

ax = fig.add_subplot(1, 2, 1, projection='3d')
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.rainbow)
fig.colorbar(surf, shrink=0.4, aspect=6)

# 使用中文
zhfont1 = matplotlib.font_manager.FontProperties \
    (fname="/home/xt/PycharmProjects/python_chinese/SourceHanSansSC-Bold.otf")

# 加标题和轴标签
ax.set_xlabel('空间轴', fontproperties=zhfont1)
ax.set_ylabel('时间轴', fontproperties=zhfont1)
ax.set_zlabel('Z ', fontproperties=zhfont1)
ax.set_title("lax-w2数值解", fontproperties=zhfont1, fontsize=20)

# 画真解与数值解误差图
Z = R + C - B
ax = fig.add_subplot(1, 2, 2, projection='3d')
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.rainbow)
fig.colorbar(surf, shrink=0.4, aspect=6)

# 使用中文
zhfont1 = matplotlib.font_manager.FontProperties \
    (fname="/home/xt/PycharmProjects/python_chinese/SourceHanSansSC-Bold.otf")

# 加标题和轴标签
ax.set_xlabel('空间轴', fontproperties=zhfont1)
ax.set_ylabel('时间轴', fontproperties=zhfont1)
ax.set_zlabel('Z ', fontproperties=zhfont1)
ax.set_title("err", fontproperties=zhfont1, fontsize=20)

plt.show()
