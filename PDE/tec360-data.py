import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
xlist = np.linspace(0, 1, 50)
ylist = np.linspace(0, 2*np.pi, 40)
X, Y = np.meshgrid(xlist, ylist)
xx = X*np.cos(Y)
print(xx.shape)
yy = X*np.sin(Y)
# print(xx)

Z = np.random.random((ylist.size, xlist.size))
fig,ax=plt.subplots()
cp = ax.contourf(xx, yy, Z)
fig.colorbar(cp) # Add a colorbar to a plot
ax.set_title('Filled Contours Plot')
#ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')

fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
surf = ax.contourf(Y, X, Z, cmap=cm.rainbow)
fig.colorbar(surf)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

ax.plot_surface(xx, yy, Z, rstride=1, cstride=1, cmap=cm.rainbow)

xx = xx.transpose()
xx = xx.reshape(xx.shape[0]*xx.shape[1],)

yy = yy.transpose()
yy = yy.reshape(yy.shape[0]*yy.shape[1],)

zz = Z.transpose()
zz = zz.reshape(Z.shape[0]*Z.shape[1],)

data = np.array([xx, yy, zz])
# print(data.shape)

data = data.transpose()

# np.savetxt('/home/xt/github/python3/file_and_animation/cmp.plt', np.c_[data],
# fmt='%.8f',delimiter='\t')



plt.show()