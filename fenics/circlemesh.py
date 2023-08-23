import dolfin
import mshr
import matplotlib.pyplot as plt

domain = mshr.Circle(dolfin.Point(0.0, 0.0), 1)
mesh = mshr.generate_mesh(domain, 10)
dolfin.plot(mesh)
plt.show()
