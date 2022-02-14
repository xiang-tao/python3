import numpy as np
from fealpy.mesh import MeshFactory as MF
import matplotlib.pyplot as plt

box2d = [0, 1, 0, 1]  # 指定x与y的最大最小值

# mesh = MF.boxmesh2d(box2d, nx=10, ny=10, meshtype='tri')
# mesh = MF.boxmesh2d(box2d, nx=10, ny=10, meshtype='quad')
# mesh = MF.boxmesh2d(box2d, nx=4, ny=4, meshtype='poly')
# mesh = MF.special_boxmesh2d(box2d, n=2, meshtype='fishbone')
# mesh = MF.special_boxmesh2d(box2d, n=10, meshtype='rice')
# mesh = MF.special_boxmesh2d(box2d, n=10, meshtype='cross')
# mesh = MF.special_boxmesh2d(box2d, n=1, meshtype='nonuniform')
# mesh = MF.unitcirclemesh(0.1,meshtype='poly')
mesh = MF.triangle(box2d,0.1,meshtype='poly')
fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
plt.show()
