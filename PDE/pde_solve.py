import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from pdedata import PdeData
from fealpy.mesh import MeshFactory as mf
from fealpy.functionspace import LagrangeFiniteElementSpace
from fealpy.boundarycondition import DirichletBC

pde = PdeData()
mesh = mf.unitcirclemesh(0.04, meshtype='tri')
space = LagrangeFiniteElementSpace(mesh, p=2)

uh = space.function()  # 返回一个有限元函数，初始自由度值是 0
# M = space.mass_matrix()
# A = space.stiff_matrix(c=pde.h_function)

qf = mesh.integrator(4, 'cell')
bcs, ws = qf.get_quadrature_points_and_weights()
cellmeasure = mesh.entity_measure('cell')
gphi = space.grad_basis(bcs)
pp = space.mesh.bc_to_point(bcs)
h_val = pde.h_function(pp)
A = np.einsum('i, ij, ijkl, ijml, j->jkm', ws, h_val, gphi, gphi, cellmeasure)

gdof = space.number_of_global_dofs()
cell2dof = space.cell_to_dof()
I = np.broadcast_to(cell2dof[:, :, None], shape=A.shape)
J = np.broadcast_to(cell2dof[:, None, :], shape=A.shape)
A = csr_matrix((A.flat, (I.flat, J.flat)), shape=(gdof, gdof))

# A = space.stiff_matrix(c=pde.h_function)

F = space.source_vector(pde.source)
bc = DirichletBC(space, pde.dirichlet)
A, F = bc.apply(A, F, uh)
uh[:] = spsolve(A, F)
print(max(uh))
print(min(uh))
errorl = space.integralalg.L2_error(pde.solution, uh.value)
print(errorl)

fig = plt.figure()
axes = fig.gca(projection='3d')
uh.add_plot(axes, cmap='rainbow')

bc = np.array([1/3, 1/3, 1/3])
val = uh(bc)
mesh.add_plot(plt, cellcolor=val, linewidths=0, showaxis=True, showcolorbar=True)

plt.show()
