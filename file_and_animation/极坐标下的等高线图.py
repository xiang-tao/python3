import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
# Using linspace so that the endpoint of 360 is included
actual = np.radians(np.linspace(0, 360, 20))
expected = np.arange(0, 70, 10)

r, theta = np.meshgrid(expected, actual)

# values表示取值大小，这里是随机产生的数字
# values = np.random.random((actual.size, expected.size))

values = 3*np.sin(theta)

fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
surf = ax.contourf(theta, r, values, cmap=cm.rainbow)
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()