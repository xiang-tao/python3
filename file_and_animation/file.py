import matplotlib.pyplot as plt
import matplotlib.animation as animation

x, y = [], []
x1, y1 = [], []

with open('data.dat', 'r') as reader:
    for line in reader.readlines():
        temp = line.split(' ')
        x.append(temp[0])
        temp1 = temp[1].split('\n')
        y.append(temp1[0])
del x[0]
del y[0]
xx = list(map(float, x))
yy = list(map(float, y))

with open('data1.dat', 'r') as reader:
    for line in reader.readlines():
        temp = line.split(' ')
        x1.append(temp[0])
        temp1 = temp[1].split('\n')
        y1.append(temp1[0])
del x1[0]
del y1[0]
xx1 = list(map(float, x1))
yy1 = list(map(float, y1))

fig, ax = plt.subplots()
xdata, ydata = [], []
xdata1, ydata1 = [], []

ax.set_xlim(-1.1, 1.1)
ax.set_ylim(0, 1)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_title("many animation lines")


def update(i):
    xdata.append(xx[i])
    ydata.append(yy[i])
    xdata1.append(xx1[i])
    ydata1.append(yy1[i])
    # ln.set_data(xdata, ydata)
    # 动画保存标签的方式如下，注意细节，在l,以及[0]，否则报错
    l, = ax.plot(xdata, ydata, 'r--', label="result")
    k = ax.plot(xdata1, ydata1, 'blue', label="real")[0]
    ax.legend([l, k], [l.get_label(), k.get_label()], loc=0)


anim = animation.FuncAnimation(fig, update, frames=len(xx), interval=40)

# anim.save('many_lines_animation.gif',writer='imagemagick')
plt.show()

# plt.plot(xx, yy, color='red', linewidth=2, linestyle='--', label='Solve the image')
# plt.plot(xx1, yy1, color='#2ca02c', linewidth=2, marker='o', label='Original image')
# plt.legend()
# plt.show()
