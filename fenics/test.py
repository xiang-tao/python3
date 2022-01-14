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

# Create mesh and define function space
mesh = UnitSquareMesh(8, 8)
print("mesh_type :", type(mesh))
V = FunctionSpace(mesh, 'P', 1)
print("V:", type(V))

# Define boundary condition
u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)
print("u_D:", type(u_D))

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Define variational problem
u = TrialFunction(V)
print("u:", type(u))
v = TestFunction(V)
f = Constant(-6.0)
a = dot(grad(u), grad(v))*dx
print(a)
print("dx", type(dx))
print("a:", type(a))
L = f*v*dx
print("L", type(L))
print("f:", type(f))
print("v:", type(v))

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Plot solution and mesh
plot(u)
# plot(mesh)

# Hold plot
plt.show()