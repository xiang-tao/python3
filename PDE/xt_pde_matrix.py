"""
该文件主要用于生成求解偏微分方程时候对应的线性方程组矩阵，维度包括一维和二维

"""
import numpy as np


def csrmatrix(A, data, indices, indptr):
    k = 0
    nr = A.shape[0]
    nc = A.shape[1]
    for i in range(nr):
        for j in range(nc):
            if A[i][j] != 0:
                data[k] = A[i][j]
                indices.append(j)
                k += 1
        indptr.append(k)


# 椭圆型


def pde_1(n):
    """
    该函数用于生成一个求解一维偏微分方程的矩阵
    :param n:矩阵的维度
    :return:返回一个矩阵
    """
    a = np.zeros((n, n))
    a[0][0] = 1
    a[n - 1][n - 1] = 1
    for i in range(n - 2):
        for j in range(n - 2):
            if i == j:
                a[i + 1][i + 1] = 2
                a[i][j + 1] = -1
                a[i + 2][j + 1] = -1
    a[0][1] = 0
    a[n - 1][n - 2] = 0
    return a


# print(pde_1(8))

def pde_2(N):
    a = np.zeros((N ** 2, N ** 2))
    for i in range(N ** 2):
        for j in range(N ** 2):
            if i < N:
                if i == j:
                    if i == 0:
                        a[i][i] = 4
                        a[i][i + 1] = -1
                        a[i][i + N] = -1
                    elif i == N - 1:
                        a[i][i] = 4
                        a[i][i - 1] = -1
                        a[i][i + N] = -1
                    else:
                        a[i][i] = 4
                        a[i][i - 1] = -1
                        a[i][i + 1] = -1
                        a[i][i + N] = -1
                else:
                    continue
            elif i > N ** 2 - N - 1:
                if i == j:
                    if i == N ** 2 - N:
                        a[i][i] = 4
                        a[i][i + 1] = -1
                        a[i][i - N] = -1
                    elif i == N ** 2 - 1:
                        a[i][i] = 4
                        a[i][i - 1] = -1
                        a[i][i - N] = -1
                    else:
                        a[i][i] = 4
                        a[i][i - 1] = -1
                        a[i][i + 1] = -1
                        a[i][i - N] = -1
                else:
                    continue
            else:
                if i == j:
                    for k in range(N - 2):
                        if i == (k + 1) * N:
                            a[i][i] = 4
                            a[i][i + 1] = -1
                            a[i][i - N] = -1
                            a[i][i + N] = -1
                            break
                        elif i == (k + 2) * N - 1:
                            a[i][i] = 4
                            a[i][i - 1] = -1
                            a[i][i - N] = -1
                            a[i][i + N] = -1
                            break
                        else:
                            continue
                    else:
                        a[i][i] = 4
                        a[i][i - 1] = -1
                        a[i][i + 1] = -1
                        a[i][i - N] = -1
                        a[i][i + N] = -1

                else:
                    continue
    return a


def pde_2_vary(Nx, Ny, hx, hy):
    """

    :param Nx:
    :param Ny:
    :param hx:
    :param hy:
    :return:
    """
    # 这里的Nx Ny表示的是横坐标和纵坐标的内点个数，
    # 不含边界点,hx hy表示的是横纵坐标的步长
    a = np.zeros((Nx * Ny, Nx * Ny))
    for i in range(Nx * Ny):
        for j in range(Nx * Ny):
            if i < Ny:
                if i == j:
                    if i == 0:
                        a[i][i] = 2 * (hx ** 2 + hy ** 2)
                        a[i][i + 1] = -hx ** 2
                        a[i][i + Ny] = -hy ** 2
                    elif i == Ny - 1:
                        a[i][i] = 2 * (hx ** 2 + hy ** 2)
                        a[i][i - 1] = -hx ** 2
                        a[i][i + Ny] = -hy ** 2
                    else:
                        a[i][i] = 2 * (hx ** 2 + hy ** 2)
                        a[i][i - 1] = -hx ** 2
                        a[i][i + 1] = -hx ** 2
                        a[i][i + Ny] = -hy ** 2
                else:
                    continue
            elif i > Nx * Ny - Ny - 1:
                if i == j:
                    if i == Nx * Ny - Ny:
                        a[i][i] = 2 * (hx ** 2 + hy ** 2)
                        a[i][i + 1] = -hx ** 2
                        a[i][i - Ny] = -hy ** 2
                    elif i == Nx * Ny - 1:
                        a[i][i] = 2 * (hx ** 2 + hy ** 2)
                        a[i][i - 1] = -hx ** 2
                        a[i][i - Ny] = -hy ** 2
                    else:
                        a[i][i] = 2 * (hx ** 2 + hy ** 2)
                        a[i][i - 1] = -hx ** 2
                        a[i][i + 1] = -hx ** 2
                        a[i][i - Ny] = -hy ** 2
                else:
                    continue
            else:
                if i == j:
                    for k in range(Ny - 2):
                        if i == (k + 1) * Ny:
                            a[i][i] = 2 * (hx ** 2 + hy ** 2)
                            a[i][i + 1] = -hx ** 2
                            a[i][i - Ny] = -hy ** 2
                            a[i][i + Ny] = -hy ** 2
                            break
                        elif i == (k + 2) * Ny - 1:
                            a[i][i] = 2 * (hx ** 2 + hy ** 2)
                            a[i][i - 1] = -hx ** 2
                            a[i][i - Ny] = -hy ** 2
                            a[i][i + Ny] = -hy ** 2
                            break
                        else:
                            continue
                    else:
                        a[i][i] = 2 * (hx ** 2 + hy ** 2)
                        a[i][i - 1] = -hx ** 2
                        a[i][i + 1] = -hx ** 2
                        a[i][i - Ny] = -hy ** 2
                        a[i][i + Ny] = -hy ** 2

                else:
                    continue
    return a


def pde_3(Nx, Ny):
    """
    该函数用于生成第三题的边界处理矩阵，对应为当边界为两点格式时候
    :param Nx: 一行块矩阵的个数
    :param Ny: 块矩阵的维度
    :return: 返回一个矩阵
    """
    a = np.zeros((Nx * Ny, Nx * Ny))
    k = 1
    for i in range(Nx * Ny):
        if i < Ny:
            a[i][i] = 1
            a[i][Ny + i] = -1
        elif i > Ny * Nx - Ny - 1:
            a[i][i] = 1
            a[i][i - Ny] = -1
        elif i == Ny * k:
            a[i][i] = 1
            a[i][i + 1] = -1
            k += 1
        elif i == Ny * k - 1:
            a[i][i] = 1
            a[i][i - 1] = -1
        else:
            a[i][i] = 4
            a[i][i + 1] = -1
            a[i][i - 1] = -1
            a[i][i + Ny] = -1
            a[i][i - Ny] = -1
    return a


def pde_3_vary(Nx, Ny):
    """
    该函数用于生成第三题的边界处理矩阵，对应为当边界为五点格式时候，即包含虚点的情况
    :param Nx: 一行块矩阵的个数
    :param Ny: 块矩阵的维度
    :return: 返回一个矩阵
    """
    a = np.zeros((Nx * Ny, Nx * Ny))
    k = 1
    for i in range(Nx * Ny):
        if i < Ny:
            if i == 0:
                a[i][i] = 1
                a[i][i + 1] = -2
                a[i][Ny + i] = -1
            elif i == Ny - 1:
                a[i][i] = 1
                a[i][i - 1] = -1
                a[i][Ny + i] = -2
            else:
                a[i][i] = 1
                a[i][i - 1] = -1
                a[i][i + 1] = -1
                a[i][Ny + i] = -2
        elif i > Ny * Nx - Ny - 1:
            if i == Ny * Nx - Ny:
                a[i][i] = 1
                a[i][i + 1] = -2
                a[i][i - Ny] = -1
            elif i == Ny * Nx - 1:
                a[i][i] = 1
                a[i][i - 1] = -1
                a[i][i - Ny] = -2
            else:
                a[i][i] = 1
                a[i][i - 1] = -1
                a[i][i + 1] = -1
                a[i][i - Ny] = -2
        elif i == Ny * k:
            a[i][i] = 1
            a[i][i + 1] = -2
            a[i][i - Ny] = -1
            a[i][i + Ny] = -1
            k += 1
        elif i == Ny * k - 1:
            a[i][i] = 1
            a[i][i - 1] = -2
            a[i][i - Ny] = -1
            a[i][i + Ny] = -1
        else:
            a[i][i] = 4
            a[i][i + 1] = -1
            a[i][i - 1] = -1
            a[i][i + Ny] = -1
            a[i][i - Ny] = -1
    return a


# 抛物型


# 向前差分格式矩阵
def back_diff(A, n, r):
    for i in range(n):
        if i == 0:
            A[i][i] = 2 * r + 1
            A[i][i + 1] = -r
        elif i == n - 1:
            A[i][i] = 2 * r + 1
            A[i][i - 1] = -r
        else:
            A[i][i] = 2 * r + 1
            A[i][i + 1] = -r
            A[i][i - 1] = -r


# 向前差分格式矩阵
def forward_diff(A, n, r):
    for i in range(n):
        if i == 0:
            A[i][i] = 1 - 2 * r
            A[i][i + 1] = r
        elif i == n - 1:
            A[i][i] = 1 - 2 * r
            A[i][i - 1] = r
        else:
            A[i][i] = 1 - 2 * r
            A[i][i + 1] = r
            A[i][i - 1] = r


# 二维下的六点格式矩阵（Grank-Nicolson）（隐式的）
# 形如a1*u^(k+1) = a0*u^k+tau*F(下面函数中的a0 a1与之对应)
def granknicolson_2(a0, a1, nx, ny, r1, r2):
    """

    :param a0: 矩阵
    :param a1: 矩阵
    :param nx: 横坐标轴内节点数
    :param ny: 纵轴内节点数
    :param r1: 横轴网格比
    :param r2: 纵轴网格比
    :return:
    """
    N = nx * ny
    for i in range(N):
        if i < ny:
            if i == 0:
                a0[i][i] = 1 - r1 - r2
                a0[i][i + 1] = 0.5 * r2
                a0[i][i + ny] = 0.5 * r1

                a1[i][i] = 1 + r1 + r2
                a1[i][i + 1] = -0.5 * r2
                a1[i][i + ny] = -0.5 * r1
            elif i == ny - 1:
                a0[i][i] = 1 - r1 - r2
                a0[i][i - 1] = 0.5 * r2
                a0[i][i + ny] = 0.5 * r1

                a1[i][i] = 1 + r1 + r2
                a1[i][i - 1] = -0.5 * r2
                a1[i][i + ny] = -0.5 * r1
            else:
                a0[i][i] = 1 - r1 - r2
                a0[i][i + 1] = 0.5 * r2
                a0[i][i - 1] = 0.5 * r2
                a0[i][i + ny] = 0.5 * r1

                a1[i][i] = 1 + r1 + r2
                a1[i][i + 1] = -0.5 * r2
                a1[i][i - 1] = -0.5 * r2
                a1[i][i + ny] = -0.5 * r1
        elif i > N - ny - 1:
            if i == N - ny:
                a0[i][i] = 1 - r1 - r2
                a0[i][i + 1] = 0.5 * r2
                a0[i][i - ny] = 0.5 * r1

                a1[i][i] = 1 + r1 + r2
                a1[i][i + 1] = -0.5 * r2
                a1[i][i - ny] = -0.5 * r1
            elif i == N - 1:
                a0[i][i] = 1 - r1 - r2
                a0[i][i - 1] = 0.5 * r2
                a0[i][i - ny] = 0.5 * r1

                a1[i][i] = 1 + r1 + r2
                a1[i][i - 1] = -0.5 * r2
                a1[i][i - ny] = -0.5 * r1
            else:
                a0[i][i] = 1 - r1 - r2
                a0[i][i + 1] = 0.5 * r2
                a0[i][i - 1] = 0.5 * r2
                a0[i][i - ny] = 0.5 * r1

                a1[i][i] = 1 + r1 + r2
                a1[i][i + 1] = -0.5 * r2
                a1[i][i - 1] = -0.5 * r2
                a1[i][i - ny] = -0.5 * r1
        elif i % ny == 0:
            a0[i][i] = 1 - r1 - r2
            a0[i][i + 1] = 0.5 * r2
            a0[i][i + ny] = 0.5 * r1
            a0[i][i - ny] = 0.5 * r1

            a1[i][i] = 1 + r1 + r2
            a1[i][i + 1] = -0.5 * r2
            a1[i][i + ny] = -0.5 * r1
            a1[i][i - ny] = -0.5 * r1
        elif i % ny == ny - 1:
            a0[i][i] = 1 - r1 - r2
            a0[i][i - 1] = 0.5 * r2
            a0[i][i + ny] = 0.5 * r1
            a0[i][i - ny] = 0.5 * r1

            a1[i][i] = 1 + r1 + r2
            a1[i][i - 1] = -0.5 * r2
            a1[i][i + ny] = -0.5 * r1
            a1[i][i - ny] = -0.5 * r1
        else:
            a0[i][i] = 1 - r1 - r2
            a0[i][i + 1] = 0.5 * r2
            a0[i][i - 1] = 0.5 * r2
            a0[i][i + ny] = 0.5 * r1
            a0[i][i - ny] = 0.5 * r1

            a1[i][i] = 1 + r1 + r2
            a1[i][i + 1] = -0.5 * r2
            a1[i][i - 1] = -0.5 * r2
            a1[i][i + ny] = -0.5 * r1
            a1[i][i - ny] = -0.5 * r1
    return a0, a1


# 二维下的Du Fort-Frankel差分格式（显示的）
def dffrankel(a1, nx, ny, r1, r2):
    N = nx * ny
    for i in range(N):
        if i < ny:
            if i == 0:
                a1[i][i + 1] = r2
                a1[i][i + ny] = r1
            elif i == ny - 1:
                a1[i][i - 1] = r2
                a1[i][i + ny] = r1
            else:
                a1[i][i + 1] = r2
                a1[i][i - 1] = r2
                a1[i][i + ny] = r1
        elif i > N - ny - 1:
            if i == N - ny:
                a1[i][i + 1] = r2
                a1[i][i - ny] = r1
            elif i == N - 1:
                a1[i][i - 1] = r2
                a1[i][i - ny] = r1
            else:
                a1[i][i + 1] = r2
                a1[i][i - 1] = r2
                a1[i][i - ny] = r1
        elif i % ny == 0:
            a1[i][i + 1] = r2
            a1[i][i + ny] = r1
            a1[i][i - ny] = r1
        elif i % ny == ny - 1:
            a1[i][i - 1] = r2
            a1[i][i + ny] = r1
            a1[i][i - ny] = r1
        else:
            a1[i][i + 1] = r2
            a1[i][i - 1] = r2
            a1[i][i + ny] = r1
            a1[i][i - ny] = r1
    return a1


# 二维下的Richardson差分格式（显示的）
def richardson(a1, nx, ny, r1, r2):
    N = nx * ny
    for i in range(N):
        if i < ny:
            if i == 0:
                a1[i][i] = -2*(r1+r2)
                a1[i][i + 1] = r2
                a1[i][i + ny] = r1
            elif i == ny - 1:
                a1[i][i] = -2*(r1+r2)
                a1[i][i - 1] = r2
                a1[i][i + ny] = r1
            else:
                a1[i][i] = -2*(r1+r2)
                a1[i][i + 1] = r2
                a1[i][i - 1] = r2
                a1[i][i + ny] = r1
        elif i > N - ny - 1:
            if i == N - ny:
                a1[i][i] = -2*(r1+r2)
                a1[i][i + 1] = r2
                a1[i][i - ny] = r1
            elif i == N - 1:
                a1[i][i] = -2*(r1+r2)
                a1[i][i - 1] = r2
                a1[i][i - ny] = r1
            else:
                a1[i][i] = -2*(r1+r2)
                a1[i][i + 1] = r2
                a1[i][i - 1] = r2
                a1[i][i - ny] = r1
        elif i % ny == 0:
            a1[i][i] = -2*(r1+r2)
            a1[i][i + 1] = r2
            a1[i][i + ny] = r1
            a1[i][i - ny] = r1
        elif i % ny == ny - 1:
            a1[i][i] = -2*(r1+r2)
            a1[i][i - 1] = r2
            a1[i][i + ny] = r1
            a1[i][i - ny] = r1
        else:
            a1[i][i] = -2*(r1+r2)
            a1[i][i + 1] = r2
            a1[i][i - 1] = r2
            a1[i][i + ny] = r1
            a1[i][i - ny] = r1
    return a1
