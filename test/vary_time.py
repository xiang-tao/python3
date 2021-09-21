import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

list1 = []
data = np.linspace(-np.pi,np.pi, 90)
for j in range(8):
    list = []
    for i in range(len(data)):
        list.append(np.sin(data[i])+j*0.3)
    list1.append(list)


fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = plt.plot([], [], 'r--',animated=True)


def init():
    ax.set_xlim(-np.pi,np.pi)
    ax.set_ylim(-1, 4)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('vary_times')
    return ln,


def update(i):
    xdata = data
    ydata = list1[i]
    ln.set_data(xdata, ydata)
    return ln,

anim = animation.FuncAnimation(fig, update, frames=8,interval=500,
                    init_func=init,blit=True)
plt.show()
