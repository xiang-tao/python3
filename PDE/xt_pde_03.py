import numpy as np
import matplotlib.pyplot as plt
import xt_pde_matrix as matrix
import xt_pde_tools as tools

n = int(input("请输入你要划分的节点数: "))
h = 2 / (n - 1)
a = matrix.pde_3(n, n)
a1 = matrix.pde_3(n, n)
b = np.zeros(n**2)
b1 = np.zeros(n**2)
x = np.zeros(n**2)
xx = np.zeros(n**2)
for i in range(n):
    for j in range(n):
        mx = i * h
        my = j * h
        if i == 0:
            b[j] = h*np.pi*np.cos(np.pi*mx)*np.sin(np.pi*my)
        elif i == n-1:
            b[i*n + j] = h*np.pi*np.cos(np.pi*mx)*np.sin(np.pi*my)
        else:
            if j == 0:
                b[i*n] = h*np.pi*np.sin(np.pi*mx)*np.cos(np.pi*my)
            elif j == n-1:
                b[i * n + j] = h * np.pi * np.sin(np.pi * mx) * np.cos(np.pi * my)
            else:
                b[i*n + j] = h*h*(2*np.pi**2-1)*np.sin(np.pi*mx)*np.cos(np.pi*my)

for i in range(n):
    for j in range(n):
        mx = i * h
        my = j * h
        if i == 0:
            b1[i] = 2*h*np.pi*np.cos(np.pi*mx)*np.sin(np.pi*my)+h*h*(2*np.pi**2-1)\
                   *np.sin(np.pi*mx)*np.cos(np.pi*my)
        elif i == n-1:
            b1[i * n + j] = 2 * h * np.pi * np.cos(np.pi * mx) * np.sin(np.pi * my)\
                           +h*h*(2*np.pi**2-1)*np.sin(np.pi*mx)*np.cos(np.pi*my)
        else:
            if j == 0:
                b1[i * n] = 2 * h * np.pi * np.sin(np.pi * mx) * np.cos(np.pi * my)\
                           +h*h*(2*np.pi**2-1)*np.sin(np.pi*mx)*np.cos(np.pi*my)
            elif j == n - 1:
                b1[i * n + j] = 2 * h * np.pi * np.sin(np.pi * mx) * np.cos(np.pi * my) \
                               + h * h * (2 * np.pi ** 2 - 1) * np.sin(np.pi * mx) * np.cos(np.pi * my)
            else:
                b1[i*n + j] = h*h*(2*np.pi**2-1)*np.sin(np.pi*mx)*np.cos(np.pi*my)

x1, list_times, list_err_norm, err = tools.jacobi(a, x, b, 5000, 10 ** (-4))
x11, list_times1, list_err_norm1, err1 = tools.jacobi(a1, x1, b1, 5000, 10 ** (-4))
# print(x1)

# 画图

plt.plot(list_times, list_err_norm, linewidth=2,marker='o', label='err_norm_5')

plt.plot(list_times1, list_err_norm1, color='red', linewidth=2, linestyle='--', label='err_norm_2')

plt.xlabel("times")
plt.ylabel("err_norm")
plt.title("err_norm result")

plt.legend()
plt.show()
