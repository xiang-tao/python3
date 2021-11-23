import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import mpl_toolkits.mplot3d
import scipy.sparse as sp
import scipy.sparse.linalg
import scipy.linalg as la
from matplotlib import cm
import dolfin
import mshr
dolfin.parameters["reorder_dofs_serial"] = False
dolfin.parameters["allow_extrapolation"] = True
N1 = N2 = 75
mesh = dolfin.RectangleMesh(dolfin.Point(0,0),dolfin.Point(1,1),N1,N2)
showmesh = dolfin.RectangleMesh(dolfin.Point(0,0),dolfin.Point(1,1),10,10)
V = dolfin.FunctionSpace(mesh, 'Lagrange', 1)
u = dolfin.TrialFunction(V)
v = dolfin.TestFunction(V)
a = dolfin.inner(dolfin.nabla_grad(u),dolfin.nabla_grad(v))*dolfin.dx
f1 = dolfin.Constant(1.0)
L1 = f1*v*dolfin.dx
f2 = dolfin.Expression("2*(x[0]*x[0]+x[1]*x[1])",degree=1)
L2 = f2*v*dolfin.dx
u0 = dolfin.Constant(0)


def u0_boundary(x, on_boundary):
    return on_boundary


bc = dolfin.DirichletBC(V, u0, u0_boundary)
A = dolfin.assemble(a)
b = dolfin.assemble(L1)
bc.apply(A,b)
u_sol1 = dolfin.Function(V)
dolfin.solve(A, u_sol1.vector(), b)
u_sol2 = dolfin.Function(V)
dolfin.solve(a == L2, u_sol2, bc)

dolfin.plot(u_sol2)


u_mat1 = np.array(u_sol2.vector()).reshape(N1+1,N2+1)
X, Y = np.meshgrid(np.linspace(0,1,N1+1),np.linspace(0,1,N2+1))
fig, ax = plt.subplots()
print(u_mat1.shape)
print(X.shape)
cp = ax.contourf(X, Y, u_mat1)
fig.colorbar(cp)

fig1, ax1 = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax1.plot_surface(X, Y, u_mat1, cmap=cm.rainbow)

xx = X
yy = Y
data = u_mat1

xx = xx.transpose()
xx = xx.reshape(xx.shape[0] * xx.shape[1], )

yy = yy.transpose()
yy = yy.reshape(yy.shape[0] * yy.shape[1], )

zz = data.transpose()
zz = zz.reshape(data.shape[0] * data.shape[1], )

data1 = np.array([xx, yy, zz])

data1 = data1.transpose()

# np.savetxt('/home/xt/github/python3/PDE/cmp.plt', np.c_[data1],
#            fmt='%.16f', delimiter='\t')

plt.show()



