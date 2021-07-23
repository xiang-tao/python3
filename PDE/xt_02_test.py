import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri


fig = plt.figure(figsize=plt.figaspect(0.5))

#============
# First plot
#============

u = np.linspace(-1, 1, endpoint=True, num=50)
v = np.linspace(-1, 1, endpoint=True, num=20)
u, v = np.meshgrid(u, v)
u, v = u.flatten(), v.flatten()

# This is the Mobius mapping, taking a u, v pair and returning an x, y, z
# triple
x = 2*u
y = 3*v
z = 3*(u**2+v**2)

# Triangulate parameter space to determine the triangles
tri = mtri.Triangulation(u, v)

# Plot the surface.  The triangles in parameter space determine which x, y, z
# points are connected by an edge.
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.plot_trisurf(x, y, z, triangles=tri.triangles, cmap=plt.cm.Spectral)
ax.set_zlim(0, 6)
plt.title("figure:z/3 = x^2/4 + y^2/9")

plt.show()
