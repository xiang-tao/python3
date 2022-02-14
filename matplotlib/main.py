import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FuncFormatter

x = np.linspace(0.5, 3.5, 100)
y = np.sin(x)

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

ax.xaxis.set_major_locator(MultipleLocator(2.0))
ax.yaxis.set_major_locator(MultipleLocator(1.0))

ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(4))


def minor_tick(x, pos):
    if not x % 1.0:
        return ""
    return "%.2f" % x


ax.xaxis.set_minor_formatter(FuncFormatter(minor_tick))

ax.tick_params("y", which='major', length=15, width=2.0, colors='r')

ax.tick_params(which='minor', length=5, width=1.0, labelsize=10, labelcolor='b')

ax.set_xlim(0, 6)
ax.set_ylim(0, 4)

ax.plot(x, y, c=(0.25, 0.25, 1.00), lw=2, zorder=10)
# ax.plot(x, y, c=(0.25, 0.25, 1.00), lw=2, zorder=0)


ax.grid(linestyle='-', linewidth=0.5, color='r', zorder=0)
# ax.grid(linestyle='-', linewidth=0.5, color='r', zorder=10)
# ax.grid(linestyle='--', linewidth=0.5, color='0.25', zorder=0)


plt.show()
