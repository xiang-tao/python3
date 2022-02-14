import numpy as np
import Integral


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
            J, A, valx, valy = Integral.int2D33(vertices)
            k = len(valx)
            uharry = np.zeros(k)
            for i in range(k):
                for alpha in range(len(uh_local)):
                    uharry[i] += uh_local[alpha] * self.uh.basis_fun(valx[i], valy[i], alpha, vertices)
            for i in range(k):
                error += J*A[i]*(self.u(valx[i], valy[i])-uharry[i])**2
        return np.sqrt(error)


class H1Norm(object):
    def __init__(self, grad_ux, grad_uy, uh, phi):
        self.grad_ux = grad_ux
        self.grad_uy = grad_uy
        self.grad_uh = uh
        self.phi = phi
        self.N = uh.V.mesh.N
        self.P = uh.V.mesh.generate_p()
        self.T = uh.V.mesh.generate_t()
        self.Pb = uh.V.generate_pb()
        self.Tb = uh.V.generate_tb()

    def half_norm(self):
        error1 = 0.0
        error2 = 0.0
        error = 0.0
        for n in range(self.N):
            uh_local = self.grad_uh.Vector()[self.Tb[n, :], :]
            vertices = self.P[self.T[n, :], :]
            J, A, valx, valy = Integral.int2D33(vertices)
            k = len(valx)
            uharry = np.zeros([k, 2])
            dxdy = np.zeros(2)
            for i in range(k):
                for alpha in range(len(uh_local)):
                    # uharry[i] += uh_local[alpha] * self.uh.basis_fun(valx[i], valy[i], alpha, vertices)
                    dxdy += uh_local[alpha] * self.phi.basis_fun(valx[i], valy[i], alpha, vertices)
                    uharry[i][0] = dxdy[0]
                    uharry[i][1] = dxdy[1]
            for i in range(k):
                error1 += J * A[i] * (self.grad_ux(valx[i], valy[i]) - uharry[i][0]) ** 2
                error2 += J * A[i] * (self.grad_uy(valx[i], valy[i]) - uharry[i][1]) ** 2
            error += error2 + error1
        return np.sqrt(error)

