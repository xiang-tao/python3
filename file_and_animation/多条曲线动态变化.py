# importing all necessary libraries
from itertools import count
import random
import matplotlib
from matplotlib import animation
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# add random points for each line
l1 = [random.randint(-20, 4) + (points ** 1.88) / (random.randint(13, 14))
      for points in range(0, 160, 2)]
l2 = [random.randint(0, 9) + (points ** 1.9) / (random.randint(9, 11))
      for points in range(0, 160, 2)]
l3 = [random.randint(-10, 10) - (points ** 1.4) / (random.randint(9, 12))
      for points in range(0, 160, 2)]

myvar = count(0, 3)
# subplots() function you can draw
# multiple plots in one figure
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 5))

# set limit for x and y axis
axes.set_ylim(-100, 500)
axes.set_xlim(0, 250)

# style for plotting line
plt.style.use("ggplot")

# create 5 list to get store element
# after every iteration
x1, y1, y2, y3, y4 = [], [], [], [], []


def animate(i):
    x1.append(next(myvar))
    # print(x1)
    # x1.append((l1[i]))
    y1.append((l1[i]))
    y2.append((l2[i]))
    y3.append((l3[i]))

    l, = axes.plot(x1, y1, color="red", label="one")
    k, = axes.plot(x1, y2, color="gray", label="two")
    m = axes.plot(x1, y3, color="blue", label="three")[0]
    axes.legend([l, k, m], [l.get_label(), k.get_label(), m.get_label()], loc=0)


# set ani variable to call the
# function recursively
anim = FuncAnimation(fig, animate, frames=8, interval=30)

plt.show()
