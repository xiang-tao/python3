import numpy as np
import matplotlib.pyplot as plt
import frog

left, right, N = 0, np.pi, 50
mesh = frog.intervalmesh(left, right, N)
V = frog.functionspace(mesh, "Lagrange", 1)
u = frog.trialfunction(V)
v = frog.testfunction(V)
nabla_u = frog.nabla_u(u)
nabla_v = frog.nabla_v(v)
a = frog.dot(nabla_u, nabla_v)


def c(x):
    return x


A = frog.assemble_A(a, c)


mass = frog.dot(u, v)
B = frog.assemble_mass(mass, c)
A = A+B


def f(x):
    return np.sin(x)+2*x*np.cos(x)


b = frog.assemble_b(f, v)
A, b = frog.dirichlet1d(1.0, -1.0, A, b)
y = frog.solve(A, b)
frog.draw(y, V)
plt.show()

