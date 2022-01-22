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
            pass
        else:
            pass

    def generate_tb(self):
        if self.degree == 1:
            Tb = self.mesh.generate_t()
            return Tb
        elif self.degree == 2:
            pass
        else:
            pass

    def generate_boundarynodes(self, type=-1):
        N1 = self.mesh.N1
        N2 = self.mesh.N2
        boundarynodes = np.zeros((2*(N1+N2), 2), dtype=np.int64)
        if self.degree == 1:
            for cn in range(N1+1):
                boundarynodes[cn][0] = type
                all_n = cn*(N2+1)
                boundarynodes[cn][1] = all_n
            for rn in range(1, N2+1):
                boundarynodes[N1+rn][0] = type
                all_n = (N2+1)*N1+rn
                boundarynodes[N1+rn][1] = all_n
            for cn in range(N1):
                boundarynodes[2*N1+N2-cn][0] = type
                all_n = N2+cn*(N2+1)
                boundarynodes[2*N1+N2-cn][1] = all_n
            for rn in range(N2-1):
                boundarynodes[2*(N1+N2)-1-rn][0] = type
                all_n = 1+rn
                boundarynodes[2*(N1+N2)-1-rn][1] = all_n
            return boundarynodes

        elif self.degree == 2:
            pass
        else:
            pass


class BaseBasisFunction(object):
    def __init__(self):
        pass

    def __str__(self):
        print("BaseBasisFunction:所有基函数的基类，包含了基函数导数等空间，可以在此处写一些共用的方法属性等")

    def basis_fun(self, x, y, basis_number, vertices):
        pass


class TrialFunction(BaseBasisFunction):
    def __init__(self, V):
        super().__init__()
        self.V = V

    def basis_fun(self, x, y, basis_number, vertices):
        degree = self.V.degree
        if degree == 1:
            coor_x = vertices[:, 0]
            coor_y = vertices[:, 1]
            J = abs((coor_x[1]-coor_x[0])*(coor_y[2]-coor_y[0])-
                    (coor_x[2]-coor_x[0])*(coor_y[1]-coor_y[0]))
            xh = ((coor_y[2]-coor_y[0])*(x-coor_x[0])-(coor_x[2]-coor_x[0])*(y-coor_y[0]))/J
            yh = (-(coor_y[1]-coor_y[0])*(x-coor_x[0])+(coor_x[1]-coor_x[0])*(y-coor_y[0]))/J
            if basis_number == 0:
                return -xh-yh+1
            elif basis_number == 1:
                return xh
            elif basis_number == 2:
                return yh
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

    def basis_fun(self, x, y, basis_number, vertices):
        degree = self.u.V.degree
        if degree == 1:
            coor_x = vertices[:, 0]
            coor_y = vertices[:, 1]
            J = abs((coor_x[1] - coor_x[0]) * (coor_y[2] - coor_y[0]) -
                    (coor_x[2] - coor_x[0]) * (coor_y[1] - coor_y[0]))
            xh_x = (coor_y[2] - coor_y[0]) / J
            xh_y = (coor_x[0] - coor_x[2]) / J
            yh_x = (coor_y[0]-coor_y[1]) / J
            yh_y = (coor_x[1] - coor_x[0]) / J
            if basis_number == 0:
                grad_phi = np.zeros(2)
                grad_phi[0] = -xh_x-yh_x
                grad_phi[1] = -xh_y-yh_y
                return grad_phi
            elif basis_number == 1:
                grad_phi = np.zeros(2)
                grad_phi[0] = xh_x
                grad_phi[1] = xh_y
                return grad_phi
            elif basis_number == 2:
                grad_phi = np.zeros(2)
                grad_phi[0] = yh_x
                grad_phi[1] = yh_y
                return grad_phi
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

    def basis_fun(self, x, y, basis_number, vertices):
        degree = self.V.degree
        if degree == 1:
            coor_x = vertices[:, 0]
            coor_y = vertices[:, 1]
            J = abs((coor_x[1] - coor_x[0]) * (coor_y[2] - coor_y[0]) -
                    (coor_x[2] - coor_x[0]) * (coor_y[1] - coor_y[0]))
            xh = ((coor_y[2] - coor_y[0]) * (x - coor_x[0]) - (coor_x[2] - coor_x[0]) * (y - coor_y[0])) / J
            yh = (-(coor_y[1] - coor_y[0]) * (x - coor_x[0]) + (coor_x[1] - coor_x[0]) * (y - coor_y[0])) / J
            if basis_number == 0:
                return -xh - yh + 1
            elif basis_number == 1:
                return xh
            elif basis_number == 2:
                return yh
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

    def basis_fun(self, x, y, basis_number, vertices):
        degree = self.v.V.degree
        h = vertices[1] - vertices[0]
        if degree == 1:
            coor_x = vertices[:, 0]
            coor_y = vertices[:, 1]
            J = abs((coor_x[1] - coor_x[0]) * (coor_y[2] - coor_y[0]) -
                    (coor_x[2] - coor_x[0]) * (coor_y[1] - coor_y[0]))
            xh_x = (coor_y[2] - coor_y[0]) / J
            xh_y = (coor_x[0] - coor_x[2]) / J
            yh_x = (coor_y[0] - coor_y[1]) / J
            yh_y = (coor_x[1] - coor_x[0]) / J
            if basis_number == 0:
                grad_phi = np.zeros(2)
                grad_phi[0] = -xh_x - yh_x
                grad_phi[1] = -xh_y - yh_y
                return grad_phi
            elif basis_number == 1:
                grad_phi = np.zeros(2)
                grad_phi[0] = xh_x
                grad_phi[1] = xh_y
                return grad_phi
            elif basis_number == 2:
                grad_phi = np.zeros(2)
                grad_phi[0] = yh_x
                grad_phi[1] = yh_y
                return grad_phi
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