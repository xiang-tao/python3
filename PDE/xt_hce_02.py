import numpy as np
import xt_pde_matrix as matrix
import xt_pde_tools as tools
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib
from matplotlib.ticker import LinearLocator

"""
向前差分格式求解一维抛物方程
"""
e = np.e
pi = np.pi
n = int(input("请输入x轴方向的想要划分的内节点数: "))
t = int(input("请输入时间轴方向的想要划分的内部层数: "))
h = 1.0 / (n + 1)
tau = 0.1 / (t + 1)
a = np.zeros((n, n))
B = np.zeros((t + 2, n))  # 保存数据方便画图，保留的是数值解
u0 = np.zeros(n)  # 将保存第0层的初始数据，此后迭代更新至最新获取的当前层数据
u1 = np.zeros(n)  # 获取最新数据，并传递给u0,方便继续迭代下去
f = np.zeros(n)
r = tau / (h ** 2)  # 网格比

# 下面这三个数据是因为我程序原因使用了csr压缩矩阵时候存储数据信息的向量
data = np.zeros(3 * n - 2)  # 3n-2代表矩阵中非零元素的个数
indices = []
indptr = [0]
# 下面开始赋初值
for i in range(n):
    hx = (i + 1) * h
    u0[i] = np.sin(2 * pi * hx)
    B[0][i] = u0[i]  # 将第一层的数据保存到B

matrix.forward_diff(a, n, r)
matrix.csrmatrix(a, data, indices, indptr)

for i in range(t + 1):
    ht = (i + 1) * tau
    for j in range(n):
        hx = (j + 1) * h
        f[j] = (10 + 4 * pi ** 2) * np.sin(2 * pi * hx) * e ** (10 * ht)
    b = tools.csrAx(data, indices, indptr, u0)
    u1 = b + tau * f
    u0 = u1
    B[i + 1] = u0  # 将u0赋值给B的i+1行，(这样操作在python里面认为是合法的)

C = np.zeros((t + 2, n))  # 获取真解数据
for i in range(t + 2):
    mt = i * tau
    for j in range(n):
        mx = (j + 1) * h
        C[i][j] = np.sin(2 * pi * mx) * e ** (10 * mt)

# 求最大规模范数误差E
E = tools.abs_max_mat(C - B)
print("最大规模误差为: %f" % E)

#============
# First plot
#============
# 画数值解图
# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
fig = plt.figure(figsize=plt.figaspect(0.5))
# Make data.
X = np.arange(h, 1, h)
Y = np.arange(0, 0.1 + tau, tau)
X, Y = np.meshgrid(X, Y)
R = X ** 0 + Y ** 0 - 2
Z = R + B

ax = fig.add_subplot(1, 2, 1, projection='3d')
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.rainbow)
ax.set_zlim(-3.5, 3.5)
ax.set_ylim(0, 0.12)
ax.yaxis.set_major_locator(LinearLocator(3))
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.4, aspect=6)

# 显示图形对应的等高线图
ax.contourf(X, Y, Z, zdir='z', offset=-3.5, cmap=cm.rainbow)

# 使用中文
zhfont1 = matplotlib.font_manager.FontProperties \
    (fname="/home/xt/PycharmProjects/python_chinese/SourceHanSansSC-Bold.otf")

# 加标题和轴标签
ax.set_xlabel('X 空间', fontproperties=zhfont1)
ax.set_ylabel('Y 时间', fontproperties=zhfont1)
ax.set_zlabel('Z 温度', fontproperties=zhfont1)
ax.set_title("Title 数值解结果 h=1/20 tau=1/1000", fontproperties=zhfont1, fontsize=20)
# fig.suptitle('test title', fontsize=20) # 这样也可以加标题，但当图片放大后标题位置不太美观

#============
# Second plot
#============
# 画真解与数值解误差图
Z = R + C - B
ax = fig.add_subplot(1, 2, 2, projection='3d')
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.rainbow)
ax.set_zlim(-0.05, 0.05)
ax.set_ylim(0, 0.12)
ax.zaxis.set_major_locator(LinearLocator(3))
ax.yaxis.set_major_locator(LinearLocator(3))

# shrink=0.4, aspect=5>>收缩= 0.4，长宽比= 5
fig.colorbar(surf, shrink=0.4, aspect=6)
ax.set_xlabel('X 空间', fontproperties=zhfont1)
ax.set_ylabel('Y 时间', fontproperties=zhfont1)
ax.set_zlabel('Z 温度', fontproperties=zhfont1)
ax.set_title("Title Error", fontproperties=zhfont1, fontsize=20)

plt.show()
