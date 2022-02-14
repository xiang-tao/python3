"""
FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
Test problem is chosen to give an exact solution at all nodes of the mesh.

  -Laplace(u) = f    in the unit square
            u = u_D  on the boundary

  u_D = 1 + x^2 + 2y^2
    f = -6
"""

from __future__ import print_function
from dolfin import *
import matplotlib.pyplot as plt

mesh = RectangleMesh(Point(0, 0), Point(1, 1), 1, 1)
V = FunctionSpace(mesh, 'Lagrange', 1)
u = TrialFunction(V)
v = TestFunction(V)
a = inner(nabla_grad(u), nabla_grad(v)) * dx
A = assemble(a)
f = Expression('x[0]', degree=1)
L = f*v*dx
b = assemble(L)

