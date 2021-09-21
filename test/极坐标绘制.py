import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import matplotlib


# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Create the mesh in polar coordinates and compute corresponding Z.
r = np.linspace(0, 2, 50)
p = np.linspace(0, 2*np.pi, 50)
R, P = np.meshgrid(r, p)
Z = R**0-P**0+np.sin(P)

# Express the mesh in the cartesian system.
X, Y = R*np.cos(P), R*np.sin(P)

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.rainbow)
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.contourf(X, Y, Z, zdir='z', offset=-1.3, cmap=cm.rainbow)
# Tweak the limits and add latex math labels.
ax.set_zlim(-1.5, 1.5)
ax.set_xlabel(r'$\phi_\mathrm{real}$')
ax.set_ylabel(r'$\phi_\mathrm{im}$')
ax.set_zlabel(r'$V(\phi)$')

plt.show()