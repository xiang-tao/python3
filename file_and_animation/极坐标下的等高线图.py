import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
# Using linspace so that the endpoint of 360 is included
# actual = np.radians(np.linspace(0, 360, 20)) # radians() 方法将角度转换为弧度
actual = np.linspace(0, 2*np.pi, 20)
expected = np.linspace(0, 1, 10)

r, theta = np.meshgrid(expected, actual)
# values表示取值大小，这里是随机产生的数字
values = np.random.random((actual.size, expected.size))
print(values)
# values = np.ones((actual.size, expected.size))
# values = 3*np.sin(theta)

fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
surf = ax.contourf(theta, r, values, cmap=cm.rainbow)
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()