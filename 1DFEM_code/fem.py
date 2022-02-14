import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import scipy.sparse.linalg
import Mesh as mesh
import LagrangeSpace as space
import NormError as NE
import Assemble as am
import BoundaryConditions as bc


def IntervalMesh(left, right, N):
    return mesh.IntervalMesh(left, right, N)


def FunctionSpace(mesh, elements_type="Lagrange", degree=2):
    if elements_type == "Lagrange":
        return space.LagrangeSpace(mesh, degree)
    elif elements_type == "CR":
        pass
    else:
        pass


def TrialFunction(V):
    return space.TrialFunction(V)


def TestFunction(V):
    return space.TestFunction(V)


def nabla_u(u):
    return space.NablaTrialFunction(u)


def nabla_v(v):
    return space.NablaTestFunction(v)


def assemble_matrix(V, au, av, *c):
    assemble = am.Am(V)
    if len(c) == 0:
        mark = 0
        return assemble.assemble_matrix(au, av, None, mark)
    elif len(c) == 1:
        if callable(c[0]):
            f = c[0]
            mark = 1
            return assemble.assemble_matrix(au, av, f, mark)
        elif type(c[0]) is np.ndarray:
            mark = 2
            arr = c[0]
            return assemble.assemble_matrix(au, av, arr, mark)
        else:
            return "waring wrong"
    else:
        return "waring too many"


def assemble_vector(V, f, v):
    b = am.Am(V)
    if callable(f):
        mark = 1
        return b.assemble_vector(f, v, mark)
    elif type(f) is np.ndarray:
        mark = 2
        vector_b = b.assemble_vector(f, v, mark)
        return vector_b
    else:
        return "waring cin wrong"


def Dirichlet(left, right, A, b):
    A, b = bc.Dirichlet(left, right, A, b)
    return A, b


def Left_Neumann(left, right, A, b, c, val):
    A, b = bc.Left_Neumann(left, right, A, b, c, val)
    return A, b


def Right_Neumann(left, right, A, b, c, val):
    A, b = bc.Right_Neumann(left, right, A, b, c, val)
    return A, b


def Right_Robbin(left, pb, qb, A, b, c, val):
    A, b = bc.Right_Robbin(left, pb, qb, A, b, c, val)
    return A, b


def Function(V):
    u = space.Function(V)
    return u


def solve(A, u, b):
    A = sp.csr_matrix(A)
    val = sp.linalg.spsolve(A, b)
    u.generate_value(val)
    # return u


def draw(u):
    x = u.V.generate_pb()
    y = u.Vector()
    plt.plot(x, y)


def L2Error(u, uh):
    L2 = NE.L2Norm(u, uh)
    error = L2.norm()
    return error


def H1half(grad_u, uh, phi):
    H1 = NE.H1Norm(grad_u, uh, phi)
    half_norm_error = H1.half_norm()
    return half_norm_error

