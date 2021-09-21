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

# 这里num是一种策略，方便画动态图
num = []
for i in range(len(xx)):
    num.append(i)

fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = ax.plot([], [], 'r--', animated=True)


# ln, = plt.plot([], [], 'r--', animated=True),这里ax.plot和plt.plot经测试效果一样

def init():
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(0, 1)
    return ln,


def update(i):
    xdata.append(xx[i])
    ydata.append(yy[i])
    ln.set_data(xdata, ydata)
    return ln,


anim = animation.FuncAnimation(fig, update, frames=num, interval=40,
                               init_func=init, blit=True)

# anim.save('test_animation.gif',writer='imagemagick')
plt.show()

# plt.plot(xx, yy, color='red', linewidth=2, linestyle='--', label='Solve the image')
# plt.plot(xx1, yy1, color='#2ca02c', linewidth=2, marker='o', label='Original image')
# plt.legend()
# plt.show()
