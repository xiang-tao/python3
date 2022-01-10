import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import scipy.sparse.linalg
import Mesh, Assemble, Space, Integral, BoundaryProcess


def intervalmesh(left, right, N):
    return Mesh.IntervalMesh(left, right, N)


def functionspace(mesh, elements_type, degree):
    return Space.FunctionSpace(mesh, elements_type, degree)


def trialfunction(V):
    return Space.TrialFunction(V)


def testfunction(V):
    return Space.TestFunction(V)


def nabla_u(u):
    return Space.NablaTrialFunction(u)


def nabla_v(v):
    return Space.NablaTestFunction(v)


def dot(au, av):
    return Integral.Int(au, av)


def assemble_A(a, *fun):
    def c(x):
        return 1


    if len(fun) == 1:
        c = fun[0]
    if len(fun) >= 2:
        print("waring:组装矩阵A的参数过多")
    V = a.au.u.V
    A = Assemble.Am(V)
    matrix_A = A.assemble_A(a, c)
    return matrix_A


def assemble_mass(a, *fun):
    def c(x):
        return 1


    if len(fun) == 1:
        c = fun[0]
    if len(fun) >= 2:
        print("waring:组装矩阵A的参数过多")
    V = a.au.V
    A = Assemble.Am(V)
    matrix_A = A.assemble_A(a, c)
    return matrix_A


def assemble_b(f, v):
    V = v.V
    b = Assemble.Am(V)
    vector_b = b.assemble_b(f, v)
    return vector_b


def dirichlet1d(left, right, A, b):
    A, b = BoundaryProcess.dirichlet1d(left, right, A, b)
    return A, b


def function(V):
    u = Space.Function(V)
    u = u.generate_function()
    return u


def solve(A, b):
    A = sp.csr_matrix(A)
    result = sp.linalg.spsolve(A, b)
    return result


def draw(y, V):
    x = V.generate_pb()
    plt.plot(x, y)

