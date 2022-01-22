import numpy as np


class IntMatirx(object):
    def __init__(self, au, av):
        self.au = au
        self.av = av

    def int1(self, vertices, alpha, beta):
        # end = 1
        # h = vertices[end] - vertices[0]
        # u0 = self.au.basis_fun(vertices[0], alpha, vertices)
        # v0 = self.av.basis_fun(vertices[0], beta, vertices)
        # val1 = u0 * v0
        #
        # u1 = self.au.basis_fun(vertices[end], alpha, vertices)
        # v1 = self.av.basis_fun(vertices[end], beta, vertices)
        # val2 = u1 * v1
        # return (val1 + val2) * h / 2
        end = 1
        a = vertices[0]
        b = vertices[end]
        x1 = -0.8611363116
        x2 = -0.3399810436
        x3 = -x2
        x4 = -x1
        A1 = A4 = 0.3478548451
        A2 = A3 = 0.6521451549
        val1 = (a + b) / 2 + (b - a) / 2 * x1
        val2 = (a + b) / 2 + (b - a) / 2 * x2
        val3 = (a + b) / 2 + (b - a) / 2 * x3
        val4 = (a + b) / 2 + (b - a) / 2 * x4
        u1 = self.au.basis_fun(val1, alpha, vertices)
        u2 = self.au.basis_fun(val2, alpha, vertices)
        u3 = self.au.basis_fun(val3, alpha, vertices)
        u4 = self.au.basis_fun(val4, alpha, vertices)
        v1 = self.av.basis_fun(val1, beta, vertices)
        v2 = self.av.basis_fun(val2, beta, vertices)
        v3 = self.av.basis_fun(val3, beta, vertices)
        v4 = self.av.basis_fun(val4, beta, vertices)
        val = (b - a) / 2 * (
                    A1 * u1 * v1 + A2 * u2 * v2 + A3 * u3 * v3 + A4 * u4 * v4)
        return val

    def int2(self, vertices, alpha, beta, f):
        end = 1
        a = vertices[0]
        b = vertices[end]
        x1 = -0.8611363116
        x2 = -0.3399810436
        x3 = -x2
        x4 = -x1
        A1 = A4 = 0.3478548451
        A2 = A3 = 0.6521451549
        val1 = (a+b)/2+(b-a)/2*x1
        val2 = (a + b) / 2 + (b - a) / 2 * x2
        val3 = (a + b) / 2 + (b - a) / 2 * x3
        val4 = (a + b) / 2 + (b - a) / 2 * x4
        u1 = self.au.basis_fun(val1, alpha, vertices)
        u2 = self.au.basis_fun(val2, alpha, vertices)
        u3 = self.au.basis_fun(val3, alpha, vertices)
        u4 = self.au.basis_fun(val4, alpha, vertices)
        v1 = self.av.basis_fun(val1, beta, vertices)
        v2 = self.av.basis_fun(val2, beta, vertices)
        v3 = self.av.basis_fun(val3, beta, vertices)
        v4 = self.av.basis_fun(val4, beta, vertices)
        val = (b-a)/2*(A1*u1*v1*f(val1)+A2*u2*v2*f(val2)+A3*u3*v3*f(val3)+A4*u4*v4*f(val4))
        return val
        # h = vertices[end] - vertices[0]
        # u0 = self.au.basis_fun(vertices[0], alpha, vertices)
        # v0 = self.av.basis_fun(vertices[0], beta, vertices)
        # val1 = u0 * v0 * f(vertices[0])
        #
        # u1 = self.au.basis_fun(vertices[end], alpha, vertices)
        # v1 = self.av.basis_fun(vertices[end], beta, vertices)
        # val2 = u1 * v1 * f(vertices[end])
        # return (val1 + val2) * h / 2

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