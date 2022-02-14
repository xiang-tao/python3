import numpy as np
import Integral


class Am(object):
    def __init__(self, V):
        self.V = V

    def assemble_matrix(self, au, av, c, mark):
        N = self.V.mesh.N1*2*self.V.mesh.N2
        P = self.V.mesh.generate_p()
        T = self.V.mesh.generate_t()
        Pb = self.V.generate_pb()
        Tb = self.V.generate_tb()
        A = np.zeros((Pb.shape[0], Pb.shape[0]))
        N_trial_lb = Tb.shape[1]
        N_test_lb = Tb.shape[1]
        a = Integral.IntMatirx(au, av)
        if mark == 0:
            for n in range(N):
                vertices = P[T[n, :], :]
                for alpha in range(N_trial_lb):
                    for beta in range(N_test_lb):
                        r = a.int1(vertices, alpha, beta)
                        i = Tb[n][beta]
                        j = Tb[n][alpha]
                        A[i][j] += r
            return A
        if mark == 1:
            for n in range(N):
                vertices = P[T[n, :], :]
                for alpha in range(N_trial_lb):
                    for beta in range(N_test_lb):
                        r = a.int2(vertices, alpha, beta, c)
                        i = Tb[n][beta]
                        j = Tb[n][alpha]
                        A[i][j] += r
            return A
        if mark == 2:
            for n in range(N):
                vertices = P[T[n, :], :]
                for alpha in range(N_trial_lb):
                    for beta in range(N_test_lb):
                        arrnumber = Tb[n]
                        r = a.int3(vertices, alpha, beta, c, arrnumber)
                        i = Tb[n][beta]
                        j = Tb[n][alpha]
                        A[i][j] += r
            return A

    def assemble_vector(self, f, v, mark):
        if mark == 1:
            N = self.V.mesh.N1*2*self.V.mesh.N2
            P = self.V.mesh.generate_p()
            T = self.V.mesh.generate_t()
            Pb = self.V.generate_pb()
            Tb = self.V.generate_tb()
            N_test_lb = Tb.shape[1]
            b = np.zeros(Pb.shape[0])
            for n in range(N):
                vertices = P[T[n, :], :]
                J, A, valx, valy = Integral.int2D33(vertices)
                k = len(valx)
                fun = f(valx, valy)
                for j in range(k):
                    for beta in range(N_test_lb):
                        varry = v.basis_fun(valx[j], valy[j], beta, vertices)
                        r = J*A[j]*varry*fun[j]
                        i = Tb[n][beta]
                        b[i] += r
            return b
        if mark == 2:
            N = self.V.mesh.N
            P = self.V.mesh.generate_p()
            T = self.V.mesh.generate_t()
            Pb = self.V.generate_pb()
            Tb = self.V.generate_tb()
            N_test_lb = Tb.shape[1]
            bb = np.zeros(Pb.shape[0])
            for n in range(N):
                vertices = P[T[n, :], :]
                arrnumber = Tb[n]
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
                for beta in range(N_test_lb):
                    v1 = v.basis_fun(val1, beta, vertices)
                    v2 = v.basis_fun(val2, beta, vertices)
                    v3 = v.basis_fun(val3, beta, vertices)
                    v4 = v.basis_fun(val4, beta, vertices)
                    r = (b - a) / 2 * (A1 * v1 * f[arrnumber[0]] + A2 * v2 * f[arrnumber[0]]
                                       + A3 * v3 * f[arrnumber[1]] + A4 * v4 * f[arrnumber[1]])
                    i = Tb[n][beta]
                    bb[i] += r
            return bb


