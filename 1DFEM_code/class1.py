import numpy as np
import matplotlib.pyplot as plt
import fem

left, right, N = 0, 1, 40
mesh = fem.IntervalMesh(left, right, N)
V = fem.FunctionSpace(mesh)
u = fem.TrialFunction(V)
v = fem.TestFunction(V)


def f(x):
    return x**3-4*x


def c(x):
    return x


A = fem.assemble_matrix(V, fem.nabla_u(u), fem.nabla_v(v), c)
B = fem.assemble_matrix(V, u, v, c)
A += B
b = fem.assemble_vector(V, f, v)
A, b = fem.Dirichlet(0, 1, A, b)
u = fem.Function(V)
fem.solve(A, u, b)
fem.draw(u)
plt.show()