import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from cmp_FEM_duiliu_data import CmpData
from fealpy.mesh import MeshFactory as mf
from fealpy.functionspace import LagrangeFiniteElementSpace
from fealpy.boundarycondition import DirichletBC

pde = CmpData()
mesh = mf.unitcirclemesh(0.04, meshtype='tri')
space = LagrangeFiniteElementSpace(mesh, p=2)
# fig = plt.figure()
# axes = fig.gca()
# mesh.add_plot(axes)

uh = space.function()  # 返回一个有限元函数，初始自由度值是 0
# print(type(uh))
# M = space.mass_matrix()
A = space.stiff_matrix(c=pde.h_function)

F = space.source_vector(pde.source)
bc = DirichletBC(space, pde.dirichlet)
A, F = bc.apply(A, F, uh)
uh[:] = spsolve(A, F)
print(max(uh))
print(min(uh))
# print(uh.shape)
fig = plt.figure()
axes = fig.gca(projection='3d')
uh.add_plot(axes, cmap='rainbow')
# axes.set_xlabel('X')

bc = np.array([1/3, 1/3, 1/3])
val = uh(bc)
mesh.add_plot(plt, cellcolor=val, linewidths=0, showaxis=True, showcolorbar=True)
plt.xlabel('X')

plt.show()
