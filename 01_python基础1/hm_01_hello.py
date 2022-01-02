import matplotlib.pyplot as plt
import numpy
x=numpy.arange(0.0, numpy.pi,0.01)
y1=numpy.cos(5*x)
y2=numpy.sin(5*x)
y=2*y1+y2
plt.plot(x, y)
plt.title("y=2cos(5x)+sin(5x)")
plt.show()
