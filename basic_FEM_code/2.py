import numpy as np
import scipy.sparse as sp
from scipy import integrate
import scipy.sparse.linalg


class IntervalMesh(object):
    def __init__(self):
        pass


def generate_boundarynodes():
    pass


def FE_basis_local_fun_1D(x, vertices, basis_type, basis_index, basis_der_x):
    h = vertices[1] - vertices[0]
    if basis_type == 101:
        if basis_index == 0:
            if basis_der_x == 0:
                return (vertices[1]-x)/h
            elif basis_der_x == 1:
                return -1/h
            elif basis_der_x >= 2:
                return 0
            else:
                print("warning wrong")
        elif basis_index == 1:
            if basis_der_x == 0:
                return (x-vertices[0])/h
            elif basis_der_x == 1:
                return 1/h
            elif basis_der_x >= 2:
                return 0
            else:
                print("warning wrong")
        else:
            print("warning wrong")
    elif basis_type ==102:
        pass
    else:
        pass


def Gauss_quad(Pb, n, c_fun, fun = FE_basis_local_fun_1D,
               vertices, basis_type, basis_index, basis_der_x):
    int_value = 0
    c_fun(Pb[n])*fun(Pb[n], vertices, basis_type, basis_index, basis_der_x)
    return int_value



def assemble_matrix(matrix_size, number_of_elements, number_of_local_basis_fun_trial,
                    number_of_local_basis_fun_test, P, T, Tb, Pb, c=1):
    A = np.zeros((matrix_size[0], matrix_size[1]))
    for n in range(number_of_elements):
        vertices = P[:,T[:,n]]
        for alpha in range(number_of_local_basis_fun_trial):
            for beta in range(number_of_local_basis_fun_test):
                int_value = Gauss_quad(Pb, n, c)
                i = Tb[beta][n]
                j = Tb[alpha][n]
                A[i][j] += int_value


