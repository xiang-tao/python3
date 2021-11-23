from fealpy.mesh import MeshFactory as MF
import matplotlib.pyplot as plt

print(MF.__file__)
box = [0,1,0,1]
mesh = MF.unitcirclemesh(0.1, meshtype='tri')
fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
plt.show()