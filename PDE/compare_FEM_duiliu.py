import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from cmp_FEM_duiliu_data import CmpData as CmpData1
from cmpdata import CmpData as CmpData2
from fealpy.mesh import MeshFactory as mf
from fealpy.functionspace import LagrangeFiniteElementSpace
from fealpy.boundarycondition import DirichletBC

pde1 = CmpData1()
mesh1 = mf.unitcirclemesh(0.05, meshtype='tri')
space1 = LagrangeFiniteElementSpace(mesh1, p=2)

uh1 = space1.function()  # 返回一个有限元函数，初始自由度值是 0
A1 = space1.stiff_matrix(c=pde1.h_function)
F1 = space1.source_vector(pde1.source)
bc1 = DirichletBC(space1, pde1.dirichlet)
A1, F1 = bc1.apply(A1, F1, uh1)
uh1[:] = spsolve(A1, F1)

pde2 = CmpData2()
mesh2 = mf.unitcirclemesh(0.05, meshtype='tri')
space2 = LagrangeFiniteElementSpace(mesh2, p=2)

uh2 = space2.function()  # 返回一个有限元函数，初始自由度值是 0
A2 = space2.stiff_matrix(c=pde2.h_function)
F2 = space2.source_vector(pde2.source)
bc2 = DirichletBC(space2, pde2.dirichlet)
A2, F2 = bc1.apply(A2, F2, uh2)
uh2[:] = spsolve(A2, F2)

uh = space2.function()
uh[:] = uh2-uh1

print(max(uh))
print(min(uh))

fig = plt.figure()
axes = fig.gca(projection='3d')
uh.add_plot(axes, cmap='rainbow')

bc = np.array([1/3, 1/3, 1/3])
val = uh(bc)
mesh2.add_plot(plt, cellcolor=val, linewidths=0, showaxis=True, showcolorbar=True)

plt.show()
