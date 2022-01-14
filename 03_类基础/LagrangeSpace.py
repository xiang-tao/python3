import numpy as np


class LagrangeSpace(object):
    def __init__(self, mesh, degree):
        self.mesh = mesh
        self.degree = degree

    def generate_pb(self):
        if self.degree == 1:
            Pb = self.mesh.generate_p()
            return Pb
        elif self.degree == 2:
            Nb = 2*self.mesh.N
            Pb = np.zeros((Nb + 1, 1))
            h = (self.mesh.right - self.mesh.left) / Nb
            for i in range(Nb + 1):
                Pb[i][0] = i * h
            return Pb
        else:
            pass

    def generate_tb(self):
        if self.degree == 1:
            Tb = self.mesh.generate_t()
            return Tb
        elif self.degree == 2:
            N = self.mesh.N
            Tb = np.zeros((N, 3), dtype=np.int64)
            for i in range(0, N):
                Tb[i][0] = 2*i
                Tb[i][1] = 2*i + 1
                Tb[i][2] = 2*i + 2
            return Tb
        else:
            pass


class BaseBasisFunction(object):
    def __init__(self):
        pass

    def __str__(self):
        print("BaseBasisFunction:所有基函数的基类，包含了基函数导数等空间，可以在此处写一些共用的方法属性等")

    def basis_fun(self, x, basis_number, vertices):
        pass


class TrialFunction(BaseBasisFunction):
    def __init__(self, V):
        super().__init__()
        self.V = V

    def basis_fun(self, x, basis_number, vertices):
        degree = self.V.degree
        if degree == 1:
            h = vertices[1] - vertices[0]
            if basis_number == 0:
                return (vertices[1] - x) / h
            elif basis_number == 1:
                return (x - vertices[0]) / h
            else:
                print("warning wrong")

        elif degree == 2:
            h = vertices[1] - vertices[0]
            if basis_number == 0:
                return 2*((x-vertices[0])/h)**2-3*(x-vertices[0])/h+1
            elif basis_number == 1:
                return -4*((x-vertices[0])/h)**2+4*(x-vertices[0])/h
            elif basis_number == 2:
                return 2*((x-vertices[0])/h)**2-(x-vertices[0])/h
            else:
                print("warning wrong")

        else:
            pass


class NablaTrialFunction(BaseBasisFunction):
    def __init__(self, u):
        super().__init__()
        self.u = u

    def basis_fun(self, x, basis_number, vertices):
        degree = self.u.V.degree
        if degree == 1:
            h = vertices[1] - vertices[0]
            if basis_number == 0:
                return -1 / h
            elif basis_number == 1:
                return 1 / h
            else:
                print("warning wrong")

        elif degree == 2:
            h = vertices[1] - vertices[0]
            if basis_number == 0:
                return 2*(2*x-2*vertices[0])/h**2-3/h
            elif basis_number == 1:
                return -4 * (2*x-2*vertices[0])/h**2 + 4 / h
            elif basis_number == 2:
                return 2 * (2*x-2*vertices[0])/h**2 - 1 / h
            else:
                print("warning wrong")

        else:
            pass


class TestFunction(BaseBasisFunction):
    def __init__(self, V):
        super().__init__()
        self.V = V

    def basis_fun(self, x, basis_number, vertices):
        degree = self.V.degree
        h = vertices[1] - vertices[0]
        if degree == 1:
            if basis_number == 0:
                return (vertices[1] - x) / h
            elif basis_number == 1:
                return (x - vertices[0]) / h
            else:
                print("warning wrong")

        elif degree == 2:
            h = vertices[1] - vertices[0]
            if basis_number == 0:
                return 2 * ((x - vertices[0]) / h) ** 2 - 3 * (x - vertices[0]) / h + 1
            elif basis_number == 1:
                return -4 * ((x - vertices[0]) / h) ** 2 + 4 * (x - vertices[0]) / h
            elif basis_number == 2:
                return 2 * ((x - vertices[0]) / h) ** 2 - (x - vertices[0]) / h
            else:
                print("warning wrong")

        else:
            pass


class NablaTestFunction(BaseBasisFunction):
    def __init__(self, v):
        super().__init__()
        self.v = v

    def basis_fun(self, x, basis_number, vertices):
        degree = self.v.V.degree
        h = vertices[1] - vertices[0]
        if degree == 1:
            if basis_number == 0:
                return -1 / h
            elif basis_number == 1:
                return 1 / h
            else:
                print("warning wrong")

        elif degree == 2:
            h = vertices[1] - vertices[0]
            if basis_number == 0:
                return 2 * (2 * x - 2 * vertices[0]) / h ** 2 - 3 / h
            elif basis_number == 1:
                return -4 * (2 * x - 2 * vertices[0]) / h ** 2 + 4 / h
            elif basis_number == 2:
                return 2 * (2 * x - 2 * vertices[0]) / h ** 2 - 1 / h
            else:
                print("warning wrong")
        else:
            pass


class Function(TrialFunction):
    def __init__(self, V):
        super().__init__(V)
        self.V = V
        Pb = self.V.generate_pb()
        self.n = Pb.shape[0]
        self.__arr = np.zeros(self.n)

    # def __call__(self, x, y):
    #     return self.array[x], self.array[y]

    def __str__(self):
        return "想要输出向量数据请使用.Vector方法"

    def __iter__(self):
        self.k = 0
        return self

    def __next__(self):
        if self.k < self.n:
            self.k += 1
        else:
            raise StopIteration
        return self.Vector()[self.k-1]

    def __getitem__(self, item):
        return self.Vector()[item]

    def Vector(self):
        arr = self.__arr.reshape(self.n, 1)
        return arr

    def generate_value(self, arr):
        self.__arr = arr



