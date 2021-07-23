"""
Notes
-----
任给一个或者一组重心坐标点，计算出网格中每个单元上对应的笛卡尔坐标
并计算函数在区域上的数值积分

(1/3, 1/3, 1/3)   --->  (NC, 2)

"""
import numpy as np
from fealpy.mesh import MeshFactory as MF
import matplotlib.pyplot as plt

box2d = [0, 1, 0, 1]
mesh = MF.boxmesh2d(box2d, nx=1, ny=1, meshtype='tri')
NC = mesh.number_of_cells()  # NC=2
node = mesh.entity('node')
cell = mesh.entity('cell')  # shape = (NC,3)
v0 = node[cell[:, 2]] - node[cell[:, 1]]  # (NC,2)
v1 = node[cell[:, 0]] - node[cell[:, 2]]  # (NC,2)
v2 = node[cell[:, 1]] - node[cell[:, 0]]  # (NC,2)
S = 1 / 2 * np.cross(v1, v2)  # nv = np.cross(v2, -v1), shape=(NC,)
bc = np.array([1 / 3, 1 / 3, 1 / 3], dtype=mesh.ftype)  # 重心坐标
qf = mesh.integrator(1)  # 一致加密
bcs, ws = qf.quadpts, qf.weights  # 重心坐标积分点和权重
ps1 = mesh.bc_to_point(bcs)
ps = mesh.bc_to_point(bc)

"""

Notes：
-----
用高斯积分计算如下的一个简单函数，f(x,y)=2(x+y)
显然此时该函数在该区域上的积分值为2,一共两个单元，
故其中每个单元的积分值为1

"""


def f(x, y):
    return 2*(x+y)


val = np.zeros((1, NC))
x = np.zeros((len(ps1)))
y = np.zeros((len(ps1)))
valf = np.zeros((len(ps1)))

for i in range(NC):
    x = ps1[:, i, 0]
    y = ps1[:, i, 1]
    valf = f(x, y)
    val[0, i] = S[i] * ws.reshape(1, -1) @ valf.reshape(-1, 1)

print('各个单元上的积分值:', val)
print('整个区域上的积分值：', np.sum(val))
