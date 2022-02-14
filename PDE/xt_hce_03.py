import numpy as np
import xt_pde_matrix as matrix
import xt_pde_tools as tools
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib

"""
二维情形下的抛物方程求解，差分格式对应为Rank-Nicholson
"""
e = np.e
pi = np.pi
times = 5000
tol = 10 ** (-4)
nx = int(input("请输入x轴方向的想要划分的内节点数: "))
ny = int(input("请输入y轴方向的想要划分的内节点数: "))
N = nx * ny
t = float(input("请输入最终的时间: "))
nt = int(input("请输入时间轴方向的想要划分的内部层数: "))
hx = 1.0 / (nx + 1)
hy = 1.0 / (ny + 1)
tau = t / (nt + 1)

a0 = np.zeros((N, N))
a1 = np.zeros((N, N))
B = np.zeros((ny, nx))  # 保存最后一层数据方便画图，保留的是数值解
u0 = np.zeros(N)  # 将保存第0层的初始数据，此后迭代更新至最新获取的当前层数据
u1 = np.zeros(N)  # 获取最新数据，并传递给u0,方便继续迭代下去
f1 = np.zeros(N)
f2 = np.zeros(N)
r1 = tau / (hx ** 2)  # 网格比
r2 = tau / (hy ** 2)

a0, a1 = matrix.granknicolson_2(a0, a1, nx, ny, r1, r2)

# 下面这三个数据是因为我程序原因使用了csr压缩矩阵时候存储数据信息的向量
data0 = np.zeros(5 * nx * ny - 2 * (nx + ny))  # 代表矩阵中非零元素的个数
indices0 = []
indptr0 = [0]
matrix.csrmatrix(a0, data0, indices0, indptr0)

data1 = np.zeros(5 * nx * ny - 2 * (nx + ny))  # 代表矩阵中非零元素的个数
indices1 = []
indptr1 = [0]
matrix.csrmatrix(a1, data1, indices1, indptr1)

k = 0
for i in range(nx):
    x = (i + 1) * hx
    for j in range(ny):
        y = (j + 1) * hy
        u0[k] = np.sin(pi * x) * np.sin(pi * y)
        k += 1

for i in range(nt + 1):
    k = 0
    t1 = i * tau
    t2 = (i + 1) * tau
    for m in range(nx):
        x = (m + 1) * hx
        for n in range(ny):
            y = (n + 1) * hy
            f1[k] = e ** (-t1) * np.sin(pi * x) * np.sin(pi * y) * (2 * pi ** 2 - 1)
            f2[k] = e ** (-t2) * np.sin(pi * x) * np.sin(pi * y) * (2 * pi ** 2 - 1)
            k += 1
    u0 = tools.csrAx(data0, indices0, indptr0, u0) + tau * 0.5 * (f1 + f2)
    tools.csrjacobi(data1, indices1, indptr1, u1, u0, times, tol)
    u0 = u1

# 将一维向量u0的值保存到矩阵B中，便于画图
k = 0
for i in range(nx):
    for j in range(ny):
        B[j][i] = u1[k]
        k += 1

C = np.zeros((ny, nx))  # 获取真解数据
for i in range(nx):
    mx = (i+1)*hx
    for j in range(ny):
        my = (j + 1) * hy
        C[j][i] = e**(-t)*np.sin(pi*mx)*np.sin(pi*my)

# 求最大规模范数误差E
E = tools.abs_max_mat(C - B)
print("最大规模误差为: %f" % E)

#============
# First plot
#============
# 画数值解图
fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 2, 1, projection='3d')
# Make data.
X = np.arange(hx, 1, hx)
Y = np.arange(hy, 1, hy)
X, Y = np.meshgrid(X, Y)
R = X ** 0 + Y ** 0 - 2
Z = R + B

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.rainbow)
fig.colorbar(surf, shrink=0.4, aspect=6)
zhfont1 = matplotlib.font_manager.FontProperties \
    (fname="/home/xt/PycharmProjects/python_chinese/SourceHanSansSC-Bold.otf")

# 加标题和轴标签
ax.set_xlabel('X 轴', fontproperties=zhfont1)
ax.set_ylabel('Y 轴', fontproperties=zhfont1)
ax.set_zlabel('Z 温度', fontproperties=zhfont1)
ax.set_title("Title 数值解结果 t=0.2 hx=1/20 hy=1/20", fontproperties=zhfont1, fontsize=20)

#============
# Second plot
#============
# 画真解与数值解误差图
ax = fig.add_subplot(1, 2, 2, projection='3d')
Z = R + C - B
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.rainbow)
fig.colorbar(surf, shrink=0.4, aspect=6)
ax.set_xlabel('X 轴', fontproperties=zhfont1)
ax.set_ylabel('Y 轴', fontproperties=zhfont1)
ax.set_zlabel('Z 温度', fontproperties=zhfont1)
ax.set_title("Title Error t=0.2s", fontproperties=zhfont1, fontsize=20)

plt.show()
