import numpy as np
import xt_pde_matrix as matrix
import xt_pde_tools as tools
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib
"""
二维情形下的抛物方程求解，差分格式对应为Du Fort-Frankel
由于该迭代方式一层的初始温度条件无法满足，需要两层，因此首先用
Crank-Nicholson方式求取第二层的温度数据
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
r1 = tau / (hx ** 2)  # 网格比,注意此处是Crank-Nicholson，Du Fort-Frankel：2r1 and 2r2
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
# 获取二层温度数据
k = 0
t1 = 0
t2 = tau
for m in range(nx):
    x = (m + 1) * hx
    for n in range(ny):
        y = (n + 1) * hy
        f1[k] = e ** (-t1) * np.sin(pi * x) * np.sin(pi * y)*(2*pi**2-1)
        f2[k] = e ** (-t2) * np.sin(pi * x) * np.sin(pi * y) * (2 * pi ** 2 - 1)
        k += 1
u0 = tools.csrAx(data0, indices0, indptr0, u0) + tau * 0.5*(f1+f2)
tools.csrjacobi(data1, indices1, indptr1, u1, u0, times, tol)
###################################################################################
# 下面开始借助已知的两层温度用Du Fort-Frankel格式计算剩下的温度层
u2 = np.zeros(N)
f3 = np.zeros(N)
a2 = np.zeros((N,N))
a2 = matrix.dffrankel(a2,nx,ny,2*r1,2*r2)
data2 = np.zeros(4 * nx * ny - 2 * (nx + ny))  # 代表矩阵中非零元素的个数
indices2 = []
indptr2 = [0]
matrix.csrmatrix(a2,data2,indices2,indptr2)
for i in range(nt):
    k = 0
    t1 = i * tau
    t2 = (i + 1) * tau
    t3 = (i + 2) * tau
    for m in range(nx):
        x = (m + 1) * hx
        for n in range(ny):
            y = (n + 1) * hy
            f1[k] = e ** (-t1) * np.sin(pi * x) * np.sin(pi * y) * (2 * pi ** 2 - 1)
            f2[k] = e ** (-t2) * np.sin(pi * x) * np.sin(pi * y) * (2 * pi ** 2 - 1)
            f3[k] = e ** (-t3) * np.sin(pi * x) * np.sin(pi * y) * (2 * pi ** 2 - 1)
            k += 1
    u0 = tools.csrAx(data2,indices2,indptr2,u1) - (2*(r1+r2)-1)*u0+tau/2*(f3+2*f2+f1)
    u2 = 1.0/(1+2*r1+2*r2)*u0
    u0 = u1
    u1 = u2

# 将一维向量u0的值保存到矩阵B中，便于画图
k = 0
for i in range(nx):
    for j in range(ny):
        B[j][i] = u2[k]
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
ax.set_title("Du Fort-Frankel 数值解结果 t=0.1 tau=1/1000 hx=1/20 hy=1/20", fontproperties=zhfont1, fontsize=20)

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
ax.set_title("Du Fort-Frankel Error t=0.1s", fontproperties=zhfont1,
             fontsize=20, color='blue', fontweight='bold', loc='left', y=0.9)

plt.show()
