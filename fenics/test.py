import matplotlib.pyplot as plt
import dolfin as df
import mshr
import numpy as np


r = 1
# domain = mshr.Circle(df.Point(0, 0), r) - mshr.Rectangle(
#     df.Point(0, 0), df.Point(1, 1))
# mesh = mshr.generate_mesh(domain, 10)
domain = mshr.Circle(df.Point(0, 0), r) - mshr.Polygon([
    df.Point(1, 0),
    df.Point(1,
        np.sqrt(2) / 2),
    df.Point(np.sqrt(2) / 2,
        np.sqrt(2) / 2)
    ])
mesh = mshr.generate_mesh(domain, 10)

df.plot(mesh)

plt.show()
