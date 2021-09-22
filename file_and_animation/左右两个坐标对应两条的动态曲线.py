import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def updateline(num, data, line1, data2, line2):
    line1.set_data(data[..., :num])
    line2.set_data(data2[..., :num])

    time_text.set_text("Points: %.0f" % int(num))

    return line1, line2


# generating data of 100 elements
# each for line 1
x = np.linspace(0, 2 * np.pi, 100)
y = np.sin(x)
data = np.array([x, y])

# generating data of 100 elements
# each for line 2
x2 = np.linspace(0, 2 * np.pi, 100)
y2 = np.cos(x2)
data2 = np.array([x2, y2])

# setup the formating for moving files
# Writer = animation.writers['ffmpeg']
# Writer = Writer(fps=10, metadata=dict(artist="Me"), bitrate=-1)

fig = plt.figure()
ax = fig.add_subplot(111)
l, = ax.plot([], [], 'r-', label="Sin")
ax2 = ax.twinx()
k = ax2.plot([], [], 'b-', label="Cos")[0]

ax.legend([l, k], [l.get_label(), k.get_label()], loc=0)

ax.set_xlabel("X")

# axis 1
ax.set_ylim(-1.5, 1.5)
ax.set_xlim(0, 7)

# axis 2
ax2.set_ylim(-1.5, 1.5)
ax2.set_xlim(0, 7)

plt.title('Sin and Cos')
time_text = ax.text(0.1, 0.95, "", transform=ax.transAxes,
                    fontsize=15, color='red')

# set line_animation variable to call
# the function recursively
line_animation = animation.FuncAnimation(
    fig, updateline, frames=100, fargs=(data, l, data2, k))
plt.show()