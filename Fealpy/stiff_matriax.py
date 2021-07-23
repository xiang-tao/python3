import numpy as np
from fealpy.mesh import MeshFactory as MF
from fealpy.functionspace import LagrangeFiniteElementSpace

box = [0, 1, 0, 1]
mesh = MF.boxmesh2d(box, nx=1, ny=1, meshtype='tri')
NN = mesh.number_of_nodes()
space = LagrangeFiniteElementSpace(mesh, p=1)

gdof = space.number_of_local_dofs()
# print(gdof, NN)
qf = mesh.integrator(2, 'cell')
bcs, ws = qf.get_quadrature_points_and_weights()
cellmeausure = mesh.entity_measure('cell')
gphi = space.grad_basis(bcs)  # (NQ,NC,ldof,2)

A = np.einsum('q,qcid,qcjd,c->cij', ws, gphi, gphi, cellmeausure)  # 爱因斯坦求和
cell2dof = space.cell_to_dof()  # (NC,ldof)

phi = space.basis(bcs)  # (NQ,NC,ldof)
# (NC, ldof, ldof)
M = np.einsum('q,qci,qcj,c->cij', ws, phi, phi, cellmeausure)
print('bcs:\n', bcs)
print('ws:\n', ws)


def f(p):
    pi = np.pi
    x = p[..., 0]
    y = p[..., 1]
    return np.sin(pi * p[..., 0]) * np.sin(pi * p[..., 1])


ps = mesh.bc_to_point(bcs)
val = f(ps)
bb = np.einsum('q,qc,qci,c->ci', ws, val, phi, cellmeausure)
b = np.zeros(gdof, dtype=np.float64)
np.add.at(b, cell2dof, bb)
