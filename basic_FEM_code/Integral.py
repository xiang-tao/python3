import numpy as np


class Int(object):
    def __init__(self, au, av):
        self.au = au
        self.av = av

    def dot(self, vertices, alpha, beta, c):
        h = vertices[1]-vertices[0]
        u0 = self.au.basis_fun(vertices[0], alpha, vertices)
        v0 = self.av.basis_fun(vertices[0], beta, vertices)
        val1 = u0*v0*c(vertices[0])

        u1 = self.au.basis_fun(vertices[1], alpha, vertices)
        v1 = self.av.basis_fun(vertices[1], beta, vertices)
        val2 = u1 * v1 * c(vertices[1])
        return (val1+val2)*h/2
