import numpy as np
import matplotlib.pyplot as plt
import fem
import NormError

left, right, N = 0, 1, 32
mesh = fem.IntervalMesh(left, right, N)
V = fem.FunctionSpace(mesh, "Lagrange", 1)
u = fem.TrialFunction(V)
v = fem.TestFunction(V)


def f(x):
    # return 50*np.cos(5*x)+25*np.sin(5*x)
    return -np.exp(x)*(np.cos(x)-2*np.sin(x)-x*np.cos(x)-x*np.sin(x))


def c(x):
    return np.exp(x)


def solution(x):
    return x*np.cos(x)


def dsolution(x):
    return np.cos(x) - x*np.sin(x)


Pb = V.generate_pb()
uu = solution(Pb)

A = fem.assemble_matrix(V, fem.nabla_u(u), fem.nabla_v(v), c)
b = fem.assemble_vector(V, f, v)
A, b = fem.Dirichlet(0, np.cos(1), A, b)
# print(A)
# A, b = fem.Left_Neumann(0, 1, A, b, p, left)
# A, b = fem.Right_Neumann(0, 2, A, b, p, right)
# A, b = fem.Right_Robbin(0, 4, 2, A, b, p, right)
uh = fem.Function(V)
fem.solve(A, uh, b)
# print(max(abs(u.Vector()-uu)))
# print(u.Vector())
L2error = NormError.L2Norm(solution, uh)
H1half = NormError.H1Norm(dsolution, uh, fem.nabla_u(u))
print(L2error.norm())
print(H1half.half_norm())
fem.draw(uh)
plt.show()

