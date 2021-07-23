import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import xt_pde_matrix as matrix
import xt_pde_tools as tools

dx = int(input("请输入区域中很坐标的长度: "))
dy = int(input("请输入区域中纵坐标的长度: "))

nx = int(input("请输入你要划分x轴方向的节点数: "))
hx = dx / (nx - 1)
ny = int(input("请输入你要划分y轴方向的节点数: "))
hy = dy / (ny - 1)
a = matrix.pde_2_vary(nx-2, ny-2, hx, hy)
b = np.zeros((nx-2)*(ny-2))
x = np.zeros((nx-2)*(ny-2))
# 此处要结合题的求解区域而定(例如区域的横纵坐标是否
# 都是从0开始，这样或许会有很多好处)
k = 0
for i in range(nx - 2):
    for j in range(ny-2):
        mx = (i + 1) * hx
        my = (j + 1) * hy
        if i == 0:
            b[k] = hx**2*hy**2*np.e**mx*np.sin(np.pi*my)*(np.pi**2-1)+hy**2*np.sin(2*np.pi*my)
            k = k+1
        elif i == nx-3:
            b[k] = hx ** 2 * hy ** 2 * np.e ** mx * np.sin(np.pi * my) * (np.pi ** 2 - 1) + hy**2*np.e**2*np.sin(2 * np.pi * my)
            k = k + 1
        else:
            b[k] = hx ** 2 * hy ** 2 * np.e ** mx * np.sin(np.pi * my) * (np.pi ** 2 - 1)
            k = k + 1

x1, list_times, list_err_norm, err = tools.jacobi(a, x, b, 10000, 10 ** (-8))
ec = max(err)
print(ec)
# 求解e0
e0 = 0
for i in range((nx-2)*(ny - 2)):
    e0 = e0 + hx*hy*err[i]**2
print(e0)

# 画图
plt.plot(list_times, list_err_norm, linewidth=2, label='err_norm')
plt.xlabel("times")
plt.ylabel("err_norm")
plt.title("err_norm result")

plt.legend()
plt.show()
