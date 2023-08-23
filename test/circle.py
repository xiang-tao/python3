import numpy as np
import matplotlib.pyplot as plt


a1 = 4
r1 = 4
a2 = 8
r2 = 8
theta = np.arange(0, 2 * np.pi, 0.01)
x = a1 + r1 * np.cos(theta)
y = r1 * np.sin(theta)

x1 = a2 + r2 * np.cos(theta)
y1 = r2 * np.sin(theta)
x2 = r1 * np.cos(theta)
y2 = a1 + r1 * np.sin(theta)

x3 = r2 * np.cos(theta)
y3 = a2 + r2 * np.sin(theta)
plt.plot(x, y)
plt.plot(x1, y1)
plt.plot(x2, y2)
plt.plot(x3, y3)

plt.show()
