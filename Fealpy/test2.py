"""
Notes:
-----
仿照三角形求解四面体的重心坐标函数的梯度

"""
import numpy as np
from fealpy.mesh import MeshFactory as MF

box3d = [0, 1, 0, 1, 0, 1]
mesh = MF.boxmesh3d(box3d, nx=1, ny=1, nz=1, meshtype='tet')
NC = mesh.number_of_cells()
node = mesh.entity('node')
cell = mesh.entity('cell')
v0 = node[cell[:, 0], :] - node[cell[:, 1], :]
v1 = node[cell[:, 0], :] - node[cell[:, 2], :]
v2 = node[cell[:, 0], :] - node[cell[:, 3], :]
v3 = node[cell[:, 2], :] - node[cell[:, 1], :]
v4 = node[cell[:, 3], :] - node[cell[:, 2], :]
v5 = node[cell[:, 1], :] - node[cell[:, 3], :]
mes = np.einsum('ij,ij->i', v1, np.cross(v3, v4))  # 爱因斯坦求和
Dlambda = np.zeros((NC, 4, 3), dtype=np.float64)
Dlambda[:, 0, :] = np.cross(v3, v4) / mes.reshape(-1, 1)
Dlambda[:, 1, :] = np.cross(-v4, v2) / mes.reshape(-1, 1)
Dlambda[:, 2, :] = np.cross(v2, -v0) / mes.reshape(-1, 1)
Dlambda[:, 3, :] = np.cross(-v3, v0) / mes.reshape(-1, 1)
print("Dlambda:", Dlambda)
