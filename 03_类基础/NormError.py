import numpy as np


class L2Norm(object):
    def __init__(self, u, uh):
        self.u = u
        self.uh = uh
        self.N = uh.V.mesh.N
        self.P = uh.V.mesh.generate_p()
        self.T = uh.V.mesh.generate_t()
        self.Pb = uh.V.generate_pb()
        self.Tb = uh.V.generate_tb()

    def norm(self):
        error = 0.0
        for n in range(self.N):
            uh_local = self.uh.Vector()[self.Tb[n, :], :]
            vertices = self.P[self.T[n, :], :]
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
            uh1, uh2, uh3, uh4 = 0.0, 0.0, 0.0, 0.0
            for alpha in range(len(uh_local)):
                uh1 += uh_local[alpha] * self.uh.basis_fun(val1, alpha, vertices)
                uh2 += uh_local[alpha] * self.uh.basis_fun(val2, alpha, vertices)
                uh3 += uh_local[alpha] * self.uh.basis_fun(val3, alpha, vertices)
                uh4 += uh_local[alpha] * self.uh.basis_fun(val4, alpha, vertices)

            error += (b-a)/2*(A1 * (self.u(val1) - uh1) ** 2 + A2 * (self.u(val2) - uh2) ** 2
                     + A3 * (self.u(val3) - uh3) ** 2 + A4 * (self.u(val4) - uh4) ** 2)
        return np.sqrt(error)


class H1Norm(object):
    def __init__(self, grad_u, uh, phi):
        self.grad_u = grad_u
        self.grad_uh = uh
        self.phi = phi
        self.N = uh.V.mesh.N
        self.P = uh.V.mesh.generate_p()
        self.T = uh.V.mesh.generate_t()
        self.Pb = uh.V.generate_pb()
        self.Tb = uh.V.generate_tb()

    def half_norm(self):
        error = 0.0
        for n in range(self.N):
            uh_local = self.grad_uh.Vector()[self.Tb[n, :], :]
            vertices = self.P[self.T[n, :], :]
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
            uh1, uh2, uh3, uh4 = 0.0, 0.0, 0.0, 0.0
            for alpha in range(len(uh_local)):
                uh1 += uh_local[alpha]*self.phi.basis_fun(val1, alpha, vertices)
                uh2 += uh_local[alpha]*self.phi.basis_fun(val2, alpha, vertices)
                uh3 += uh_local[alpha]*self.phi.basis_fun(val3, alpha, vertices)
                uh4 += uh_local[alpha]*self.phi.basis_fun(val4, alpha, vertices)

            error += (b - a) / 2 * (A1 * (self.grad_u(val1) - uh1) ** 2 + A2 * (self.grad_u(val2) - uh2) ** 2
                                    + A3 * (self.grad_u(val3) - uh3) ** 2 + A4 * (self.grad_u(val4) - uh4) ** 2)
        return np.sqrt(error)
