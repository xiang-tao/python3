import numpy as np


class Am(object):
    def __init__(self, V):
        self.V = V

    def assemble_A(self, a, c):
        N = self.V.mesh.N
        P = self.V.mesh.generate_p()
        T = self.V.mesh.generate_t()
        Tb = self.V.generate_tb()
        A = np.zeros((N + 1, N + 1))
        N_trial_lb = Tb.shape[1]
        N_test_lb = Tb.shape[1]
        for n in range(N):
            vertices = P[T[n, :], :]
            verticesarry = T[n]
            for alpha in range(N_trial_lb):
                for beta in range(N_test_lb):
                    r = a.dot(vertices, alpha, beta, c, verticesarry)
                    i = Tb[n][beta]
                    j = Tb[n][alpha]
                    A[i][j] += r
        return A

    def assemble_b(self, f, v):
        N = self.V.mesh.N
        P = self.V.mesh.generate_p()
        T = self.V.mesh.generate_t()
        Tb = self.V.generate_tb()
        N_test_lb = Tb.shape[1]
        b = np.zeros(N + 1)
        for n in range(N):
            vertices = P[T[n, :], :]
            for beta in range(N_test_lb):
                h = vertices[1] - vertices[0]
                v0 = v.basis_fun(vertices[0], beta, vertices)
                val0 = f[T[n][0]]*v0
                v1 = v.basis_fun(vertices[1], beta, vertices)
                val1 = f[T[n][1]] * v1
                r = (val0+val1)*h/2
                i = Tb[n][beta]
                b[i] += r
        return b
