import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import mpl_toolkits.mplot3d
import scipy.sparse as sp
import scipy.sparse.linalg
import scipy.linalg as la
from matplotlib import cm
import dolfin
import mshr

r = 1
domain = mshr.Circle(dolfin.Point(0, 0), r)
mesh = mshr.generate_mesh(domain, 10)
# dolfin.plot(mesh)
V = dolfin.FunctionSpace(mesh, 'Lagrange', 1)


def u0_boundary(x, on_boundary):
    x, y = x[0], x[1]
    return on_boundary and abs(np.sqrt(x ** 2 + y ** 2) - r) < 5e-2


bcs = dolfin.DirichletBC(V, dolfin.Constant(10), u0_boundary)

u = dolfin.TrialFunction(V)
v = dolfin.TestFunction(V)
a = dolfin.inner(dolfin.nabla_grad(u), dolfin.nabla_grad(v)) * dolfin.dx
f = dolfin.Constant(0.0)
L = f * v * dolfin.dx
u_sol = dolfin.Function(V)

dolfin.solve(a == L, u_sol, bcs)
dolfin.plot(u_sol)

# coordinates = mesh.coordinates()
# triangles = mesh.cells()
# triangulation = mpl.tri.Triangulation(coordinates[:,0], coordinates[:,1],triangles)
# fig = (ax1,ax2) = plt.subplots(1,2)
# plt.triplot(triangulation)
# c= plt.tripcolor(triangulation,np.array(u_sol.vector()),cmap = mpl.cm.get_cmap("Reds"))
# # cb = plt.colorbar(c,ax = ax2)

plt.show()
