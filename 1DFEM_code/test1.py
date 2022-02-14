import numpy as np
import matplotlib.pyplot as plt
import fem

left, right, N = -1, 1, 64
mesh = fem.IntervalMesh(left, right, N)
V = fem.FunctionSpace(mesh, degree=2)
u = fem.TrialFunction(V)
v = fem.TestFunction(V)


def f(x):
    return np.exp(-x*x)*(4*x**4-15*x**2+5)


def solution(x):
    return np.exp(-x*x)*(1-x**2)


def dsolution(x):
    return np.exp(-x**2)*(2*x**3-4*x)


A = fem.assemble_matrix(V, fem.nabla_u(u), fem.nabla_v(v))
B = fem.assemble_matrix(V, u, v)
A += B
b = fem.assemble_vector(V, f, v)
A, b = fem.Dirichlet(0, 0, A, b)
uh = fem.Function(V)
fem.solve(A, uh, b)
L2error = fem.L2Error(solution, uh)
H1halferror = fem.H1half(dsolution, uh, fem.nabla_u(u))
print(L2error)
print(H1halferror)
fem.draw(uh)
plt.show()