import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm

list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax.set_zlim(-2, 2)

def update(i, fig):
    X = np.arange(-5, 5, 0.25)
    Y = np.arange(-5, 5, 0.25)
    X, Y = np.meshgrid(X, Y)
    R = np.sqrt(X ** 2 + Y ** 2)
    # add random shift in z directions
    Z = np.sin(R) * i * 0.2
    # Z = np.sin(R) + np.random.random_sample()
    if len(fig.axes[0].collections) != 0:
        fig.axes[0].collections = []
        surf = fig.axes[0].plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    else:
        surf = fig.axes[0].plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.set_title('vary_times' + list[i])
    fig.canvas.draw()
    time_text.set_text("Times: %.1f" % (0.5 * i))
    return surf,

time_text = ax.text(0, 0, 0, "", transform=ax.transAxes,
                     fontsize=15, color='red')
anim = animation.FuncAnimation(fig, update, fargs=[fig], frames=10, interval=500, blit=True)
# anim.save('3d_vary.gif', writer='imagemagick')
plt.show()
