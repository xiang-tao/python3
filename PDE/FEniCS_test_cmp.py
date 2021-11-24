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
sin = np.sin
cos = np.cos
pi = np.pi
dolfin.parameters["reorder_dofs_serial"] = False
dolfin.parameters["allow_extrapolation"] = True

ww = 5.0 * pi / 3.0  # 50.0  # ----------------------------------------晶片的角速度
wp = 10.0 * pi / 3.0  # 100.0  # ----------------------------------------抛光垫的角速度
angle1 = 0.02 * pi / 180  # ---------------------------转角
angle2 = 0.018 * pi / 180  # ---------------------------倾角
d = 150000.0  # ---------------------------晶片和抛光垫的旋转中心距
r0 = 50000.0
p0 = 101000.0
hpiv = 80.0  # ----------------------------------------晶片中心高度
viscosity = 0.00214  # ---------------------------抛光液粘度
k = hpiv**3/r0*p0/r0
kk = r0/hpiv

r_inner = 1
domain = mshr.Circle(dolfin.Point(0, 0), r_inner)
mesh = mshr.generate_mesh(domain, 10)

V = dolfin.FunctionSpace(mesh, 'Lagrange', 1)
u = dolfin.TrialFunction(V)
v = dolfin.TestFunction(V)

h = dolfin.Expression("1-50000.0/80.0*x[0]*sin(0.02 * pi / 180)-50000.0/80.0*x[1]"
                      "*sin(0.018 * pi / 180)", degree=1)
a = h*h*h*dolfin.inner(dolfin.nabla_grad(u),dolfin.nabla_grad(v))*dolfin.dx
# f1 = dolfin.Constant(1)
# L1 = f1*v*dolfin.dx
f2 = dolfin.Expression("6*0.00214*(((10.0 * pi / 3.0)*(x[0]*50000.0+150000.0)+"
                       "(5.0 * pi / 3.0)*50000.0*x[0])"
                       "*sin(0.02 * pi / 180)-(((5.0 * pi / 3.0)+(10.0 * pi / 3.0))"
                       "*50000.0*x[1])*sin(0.018 * pi / 180))/(80.0*80.0*80.0/"
                       "(50000.0*101000.0)/50000.0)", degree=1)
L2 = f2*v*dolfin.dx
u0 = dolfin.Constant(1)


def u0_boundary(x, on_boundary):
    return on_boundary


bc = dolfin.DirichletBC(V, u0, u0_boundary)
u_sol2 = dolfin.Function(V)
dolfin.solve(a == L2, u_sol2, bc)

dolfin.plot(u_sol2,cmap=cm.rainbow)

u_mat = np.array(u_sol2.vector())
print(max(u_mat))

plt.show()



