import numpy as np


class FunctionSpace(object):
    def __init__(self, mesh, elements_type, degree):
        self.mesh = mesh
        self.elements_type = elements_type
        self.degree = degree

    def generate_pb(self):
        if self.elements_type == "Lagrange":
            if self.degree == 1:
                Pb = self.mesh.generate_p()
                return Pb
            elif self.degree == 2:
                pass
            else:
                pass

        elif self.elements_type == "CR":
            pass

        else:
            pass

    def generate_tb(self):
        if self.elements_type == "Lagrange":
            if self.degree == 1:
                Tb = self.mesh.generate_t()
                return Tb
            elif self.degree == 2:
                pass
            else:
                pass

        elif self.elements_type == "CR":
            pass

        else:
            pass


class TrialFunction(object):
    def __init__(self, V):
        self.V = V

    def basis_fun(self, x, basis_number, vertices):
        elements_type = self.V.elements_type
        degree = self.V.degree
        h = vertices[1] - vertices[0]
        if elements_type == "Lagrange":
            if degree == 1:
                if basis_number == 0:
                    return (vertices[1] - x) / h
                elif basis_number == 1:
                    return (x - vertices[0]) / h
                else:
                    print("warning wrong")

            elif degree == 2:
                pass

            else:
                pass

        elif elements_type == "CR":
            pass

        else:
            pass


class NablaTrialFunction(object):
    def __init__(self, u):
        self.u = u

    def basis_fun(self, x, basis_number, vertices):
        elements_type = self.u.V.elements_type
        degree = self.u.V.degree
        h = vertices[1] - vertices[0]
        if elements_type == "Lagrange":
            if degree == 1:
                if basis_number == 0:
                    return -1 / h
                elif basis_number == 1:
                    return 1 / h
                else:
                    print("warning wrong")

            elif degree == 2:
                pass

            else:
                pass

        elif elements_type == "CR":
            pass

        else:
            pass


class TestFunction(object):
    def __init__(self, V):
        self.V = V

    def basis_fun(self, x, basis_number, vertices):
        elements_type = self.V.elements_type
        degree = self.V.degree
        h = vertices[1] - vertices[0]
        if elements_type == "Lagrange":
            if degree == 1:
                if basis_number == 0:
                    return (vertices[1] - x) / h
                elif basis_number == 1:
                    return (x - vertices[0]) / h
                else:
                    print("warning wrong")

            elif degree == 2:
                pass

            else:
                pass

        elif elements_type == "CR":
            pass

        else:
            pass


class NablaTestFunction(object):
    def __init__(self, v):
        self.v = v

    def basis_fun(self, x, basis_number, vertices):
        elements_type = self.v.V.elements_type
        degree = self.v.V.degree
        h = vertices[1] - vertices[0]
        if elements_type == "Lagrange":
            if degree == 1:
                if basis_number == 0:
                    return -1 / h
                elif basis_number == 1:
                    return 1 / h
                else:
                    print("warning wrong")

            elif degree == 2:
                pass

            else:
                pass

        elif elements_type == "CR":
            pass

        else:
            pass


class Function(object):
    def __init__(self, V):
        self.V = V

    def generate_function(self):
        Pb = self.V.generate_pb()
        end = Pb.shape[0]
        u = np.zeros(end + 1)
        return u
