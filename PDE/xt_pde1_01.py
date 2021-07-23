import numpy as np
import matplotlib.pyplot as plt
# import chineseize_matplotlib　
# 这是一个显示中文的包，但是出现了意料之外的问题，暂时放弃使用
import xt_pde_matrix as matrix
import xt_pde_tools as tools

n = int(input("请输入你要划分的节点数: "))
h = 2 / (n - 1)
a = matrix.pde_1(n)
data = np.zeros(3*(n-2))
indices = []
indptr = [0]
matrix.csrmatrix(a,data,indices,indptr)
b = np.zeros(n)
x = np.zeros(n)
for i in range(n - 2):
    m = ((i + 1) * h - 1)
    b[i + 1] = h * h * (4 * m ** 4 - 14 * m ** 2 + 4) * np.e ** (-m ** 2)
    # print(b[i+1])

x1, list_times, list_err_norm, err = tools.csrjacobi(data,indices,indptr, x, b, 10000, 10 ** (-8))
# print(x1)

# 画图
x2 = np.linspace(-1, 1, n)
y = np.e ** (-x2 ** 2) * (1 - x2 ** 2)

plt.subplot(211)
plt.plot(x2, y, color='#2ca02c', linewidth=2, marker='o', label='Original image')
plt.plot(x2, x1, color='red', linewidth=2, linestyle='--', label='Solve the image')
plt.xlabel("x_value", x=0.9)  # x=0.9表示标签位置
plt.ylabel("y_value")
plt.title('Compare result')
plt.legend()

plt.subplot(212)
plt.plot(list_times, list_err_norm, linewidth=2, label='err_norm')
plt.xlabel("times")
plt.ylabel("err_norm")
plt.title("err_norm result", color='blue', fontsize='15', fontweight='bold', loc='left', y=1)

plt.legend()
plt.show()
