import xt_pde_matrix as matrix
import xt_pde_tools as tools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
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
list = []
sublist = []
r = tau / (h ** 2)

hx = 0
for i in range(n):
    hx = (i + 1) * h
    u0[i] = e ** (-100 * (hx - 0.25) ** 2) + 0.1 * np.sin(20 * pi * hx)

matrix.back_diff(A, n, r)
matrix.csrmatrix(A, data, indices, indptr)
for i in range(n):
    B[0][i] = u0[i]
    sublist.append(u0[i])
list.append(sublist)
for i in range(t + 1):
    tools.csrjacobi(data, indices, indptr, u1, u0, times, tol)
    sublist = []
    for k in range(n):
        u0[k] = u1[k]
        sublist.append(u0[k])
        B[i + 1][k] = u0[k]
    list.append(sublist)
xlist = np.arange(h, 1, h)

fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = plt.plot([], [], 'r--')


def init():
    ax.set_xlim(0, 1.1)
    ax.set_ylim(0, 1)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('vary_times')

    return ln,


def update(i):
    xdata = xlist
    ydata = list[i]
    ln.set_data(xdata, ydata)
    # 有时候这个text在该函数下不起作用，
    # 那么需要将animation.FuncAnimation中的参数blit=True去除，
    # 或者把该text放到init函数下，或者放在函数外面，但是这样无法
    # 达到动态的效果，具体不起作用原因原因我也不清楚
    time_text.set_text("Times: %.3f" % (tau * i))
    return ln,


time_text = ax.text(0.4, 0.6, "", transform=ax.transAxes,
                    fontsize=15, color='red')

anim = animation.FuncAnimation(fig, update, frames=len(list), interval=50,
                               init_func=init)
plt.show()

