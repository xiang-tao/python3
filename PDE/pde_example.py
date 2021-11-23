import numpy as np
from fealpy.mesh import MeshFactory as MF
from fealpy.functionspace import LagrangeFiniteElementSpace
import matplotlib.pyplot as plt

mesh = MF.unitcirclemesh(0.5, meshtype='tri')
# mesh = MF.boxmesh2d([0,1,0,1], nx=1, ny=1, meshtype='tri')
face = mesh.number_of_faces()
cell = mesh.number_of_cells()
print(face,cell)
cellmeasure = mesh.entity_measure("cell")
print(cellmeasure)
space = LagrangeFiniteElementSpace(mesh,1)
ldof = space.number_of_local_dofs()
gdolf = space.number_of_global_dofs()
bc = np.array([1,0,0])
ps = mesh.bc_to_point(bc)
phi = space.basis(bc)
gphi = space.grad_basis(bc)
ipoints = space.interpolation_points()
print(ipoints)
print(phi)
print(gphi)
print(ps)
print(ldof,gdolf)
fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_node(axes, showindex=True)
mesh.find_edge(axes, showindex=True)
mesh.find_cell(axes, showindex=True)
plt.show()

