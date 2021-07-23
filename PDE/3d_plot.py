import xt_pde_matrix as matrix
import xt_pde_tools as tools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib
from matplotlib.ticker import LinearLocator
"""
向后差分格式求解一维抛物方程代码
"""
e = np.e
pi = np.pi
n = 19;
times = 3000;
tol = 10 ** (-4)
h = 1.0 / (n + 1)
t = 99;
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

#
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
X = np.arange(h, 1, h)
Y = np.arange(0, 0.1 + tau, tau)
X, Y = np.meshgrid(X, Y)
R = X ** 0 + Y ** 0 - 2
Z = R + B

# Plot the surface.
# cmap=cm.rainbow: 图形上的颜色分布 rainbow彩虹色，coolwarm:也是颜色种类，
# rstride=1, cstride=1:x y轴的颜色分隔宽度，即线条分隔宽度，每个线条对应一种颜色
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.rainbow)
# Customize the z axis.
ax.set_zlim(-1.01, 1.01) # 限制z轴的范围 xlim限制x轴，ylim限制y轴
# ax.zaxis.set_major_locator(LinearLocator(10)) # 将z轴平均分10段，显示10个刻度
# A StrMethodFormatter is used automatically
# ax.zaxis.set_major_formatter('{x:.02f}')

# Add a color bar which maps values to colors.
# 在图形旁边增加一个竖着的框，显示不同颜色对应的数据大小
fig.colorbar(surf, shrink=0.5, aspect=5)

# 显示图形对应的登高线图
# zidr='z'表示与z轴垂直的登高线图，即从z轴上方上往下看，同理可以有x　y轴方向的zidr='x' ,zidr='y'
# offset=-1,若是与z轴垂直的当高线，则表示在z=-1处显示登高线图
# cmap=cm.rainbow:表示彩虹色
ax.contourf(X, Y, Z, zdir='z', offset=-1, cmap=cm.rainbow)

# 使用中文
# fname 为 你下载的字体库路径，注意 SourceHanSansSC-Bold.otf 字体的路径
# 注意使用时候要在后面加　fontproperties=zhfont1，例如：ax.set_title("Title 标题",fontproperties=zhfont1)
zhfont1 = matplotlib.font_manager.FontProperties\
    (fname="/home/xt/PycharmProjects/python_chinese/SourceHanSansSC-Bold.otf")

# 加标题和轴标签，使用：fontsize=20(any number)可以改变标签字体大小
ax.set_xlabel('X 空间', fontsize=20)
ax.set_ylabel('Y 时间')
ax.set_zlabel('Z 温度')
ax.set_title("Title 标题",fontproperties=zhfont1, fontsize=20)
# fig.suptitle('test title', fontsize=20) # 这样也可以加标题，但当图片放大后标题位置不太美观

plt.show()
