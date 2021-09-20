import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from matplotlib import cm

fig = plt.figure()
ax = fig.gca(projection='3d')


def update(frame, fig):
    X = np.arange(-5, 5, 0.25)
    Y = np.arange(-5, 5, 0.25)
    X, Y = np.meshgrid(X, Y)
    R = np.sqrt(X**2 + Y**2)
    # add random shift in z directions
    Z = np.sin(R)+np.random.random_sample()
    if len(fig.axes[0].collections) != 0:
        fig.axes[0].collections = []
        surf = fig.axes[0].plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    else:
        surf = fig.axes[0].plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.set_zlim(-1.5, 1.5)

    fig.canvas.draw()
    return surf,


anim = animation.FuncAnimation(fig, update, fargs=[fig], frames=5, interval=500, blit=True)
plt.show()