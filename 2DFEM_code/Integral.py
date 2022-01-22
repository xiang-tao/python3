import numpy as np


def int2D22(vertices):
    a = np.sqrt(3)
    n = vertices.shape[0]
    x = np.zeros(n)
    y = np.zeros(n)
    valx = np.zeros(4)
    valy = np.zeros(4)
    gauss_x = np.array([(1 / a + 1) / 2, (1 / a + 1) / 2, (-1 / a + 1) / 2, (-1 / a + 1) / 2])
    gauss_y = np.array([(1 - 1 / a) * (1 + 1 / a) / 4, (1 - 1 / a) * (1 - 1 / a) / 4,
                        (1 + 1 / a) * (1 + 1 / a) / 4, (1 - 1 / a) * (1 + 1 / a) / 4])
    A = np.array([(1 - 1 / a) / 8, (1 - 1 / a) / 8, (1 + 1 / a) / 8, (1 + 1 / a) / 8])

    for i in range(n):
        x[i] = vertices[i][0]
        y[i] = vertices[i][1]

    J = abs((x[1] - x[0]) * (y[2] - y[0]) - (y[1] - y[0]) * (x[2] - x[0]))
    for i in range(4):
        valx[i] = x[0] + (x[1] - x[0]) * gauss_x[i] + (x[2] - x[0]) * gauss_y[i]
        valy[i] = y[0] + (y[1] - y[0]) * gauss_x[i] + (y[2] - y[0]) * gauss_y[i]
    return J, A, valx, valy


def int2D33(vertices):
    n = vertices.shape[0]
    x = np.zeros(n)
    y = np.zeros(n)
    valx = np.zeros(9)
    valy = np.zeros(9)
    gauss_x = np.array([0.112701665, 0.112701665, 0.112701665, 0.500000000, 0.500000000,
                        0.500000000, 0.887298334, 0.887298334, 0.887298334])
    gauss_y = np.array([0.100000000, 0.443649167, 0.787298334, 0.056350832, 0.250000000,
                        0.443649167, 0.012701665, 0.056350832, 0.100000000])
    A = np.array([0.068464377, 0.109543004, 0.068464377, 0.061728395, 0.098765432,
                  0.061728359, 0.008696116, 0.013913785, 0.008696116])
    for i in range(n):
        x[i] = vertices[i][0]
        y[i] = vertices[i][1]
    J = abs((x[1] - x[0]) * (y[2] - y[0]) - (y[1] - y[0]) * (x[2] - x[0]))
    for i in range(9):
        valx[i] = x[0] + (x[1] - x[0]) * gauss_x[i] + (x[2] - x[0]) * gauss_y[i]
        valy[i] = y[0] + (y[1] - y[0]) * gauss_x[i] + (y[2] - y[0]) * gauss_y[i]
    return J, A, valx, valy


class IntMatirx(object):
    def __init__(self, au, av):
        self.au = au
        self.av = av

    def int1(self, vertices, alpha, beta):
        J, A, valx, valy = int2D33(vertices)
        n = len(valx)
        uv = np.zeros(n)
        for i in range(n):
            u = self.au.basis_fun(valx[i], valy[i], alpha, vertices)
            v = self.av.basis_fun(valx[i], valy[i], beta, vertices)
            uv[i] = np.dot(u, v)
        r = J*np.dot(A, uv)
        return r


    def int2(self, vertices, alpha, beta, f):
        coor_x0 = vertices[0][0]
        coor_x1 = vertices[1][0]
        coor_x2 = vertices[2][0]
        coor_y0 = vertices[0][1]
        coor_y1 = vertices[1][1]
        coor_y2 = vertices[2][1]
        a = np.sqrt(3)
        x1 = (1 / a + 1) / 2
        x2 = x1
        x3 = (-1 / a + 1) / 2
        x4 = x3
        y1 = (1 - 1 / a) * (1 + 1 / a) / 4
        y2 = (1 - 1 / a) * (1 - 1 / a) / 4
        y3 = (1 + 1 / a) * (1 + 1 / a) / 4
        y4 = (1 - 1 / a) * (1 + 1 / a) / 4
        A1 = (1 - 1 / a) / 8
        A2 = (1 - 1 / a) / 8
        A3 = (1 + 1 / a) / 8
        A4 = (1 + 1 / a) / 8
        J = abs((coor_x1 - coor_x0) * (coor_y2 - coor_y0) - (coor_y1 - coor_y0) * (coor_x2 - coor_x0))
        valx1 = coor_x0 + (coor_x1 - coor_x0) * x1 + (coor_x2 - coor_x0) * y1
        valy1 = coor_y0 + (coor_y1 - coor_y0) * x1 + (coor_y2 - coor_y0) * y1
        valx2 = coor_x0 + (coor_x1 - coor_x0) * x2 + (coor_x2 - coor_x0) * y2
        valy2 = coor_y0 + (coor_y1 - coor_y0) * x2 + (coor_y2 - coor_y0) * y2
        valx3 = coor_x0 + (coor_x1 - coor_x0) * x3 + (coor_x2 - coor_x0) * y3
        valy3 = coor_y0 + (coor_y1 - coor_y0) * x3 + (coor_y2 - coor_y0) * y3
        valx4 = coor_x0 + (coor_x1 - coor_x0) * x4 + (coor_x2 - coor_x0) * y4
        valy4 = coor_y0 + (coor_y1 - coor_y0) * x4 + (coor_y2 - coor_y0) * y4
        u1 = self.au.basis_fun(valx1, valy1, alpha, vertices)
        u2 = self.au.basis_fun(valx2, valy2, alpha, vertices)
        u3 = self.au.basis_fun(valx3, valy3, alpha, vertices)
        u4 = self.au.basis_fun(valx4, valy4, alpha, vertices)
        v1 = self.av.basis_fun(valx1, valy1, beta, vertices)
        v2 = self.av.basis_fun(valx2, valy2, beta, vertices)
        v3 = self.av.basis_fun(valx3, valy3, beta, vertices)
        v4 = self.av.basis_fun(valx4, valy4, beta, vertices)
        val = J * (A1 * np.dot(u1, v1)*f(valx1, valy1) + A2 * np.dot(u2, v2)*f(valx2, valy2) + A3 * np.dot(u3, v3)*f(valx3, valy3) + A4 * np.dot(u4, v4)*f(valx4,valy4))
        return val

    def int3(self, vertices, alpha, beta, arr, arrnumber):
        end = 1
        h = vertices[end] - vertices[0]
        u0 = self.au.basis_fun(vertices[0], alpha, vertices)
        v0 = self.av.basis_fun(vertices[0], beta, vertices)
        val1 = u0 * v0 * arr[arrnumber[0]]

        u1 = self.au.basis_fun(vertices[end], alpha, vertices)
        v1 = self.av.basis_fun(vertices[end], beta, vertices)
        val2 = u1 * v1 * arr[arrnumber[end]]
        return (val1 + val2) * h / 2


class IntVector(object):
    def __init__(self, av):
        self.av = av

    def int1(self):
        pass