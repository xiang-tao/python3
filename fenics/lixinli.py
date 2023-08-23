"""
FEniCS tutorial demo program: Deflection of a membrane.

  -Laplace(w) = p  in the unit circle
            w = 0  on the boundary

The load p is a Gaussian function centered at (0, 0.6).
"""

# from __future__ import print_function
from matplotlib import cm
from fenics import *
# from dolfin import *
from mshr import *
import numpy as np

sin = np.sin
cos = np.cos
pi = np.pi
# 请注意r/s单位不是国际单位，应换算成弧度每秒，而一转的弧度是2*pi，故需要乘2*pi
ww = 5.0 * pi / 3.0  # 50.0  # ----------------------------------------晶片的角速度
wp = 10.0 * pi / 3.0  # 100.0  # ----------------------------------------抛光垫的角速度
angle1 = 0.02 * pi / 180  # ---------------------------转角
angle2 = 0.018 * pi / 180  # ---------------------------倾角
d = 0.15  # ---------------------------晶片和抛光垫的旋转中心距
r0 = 5.0 * 1e-2
p0 = 101000.0
hpiv = 8.0 * 1e-5  # ----------------------------------------晶片中心高度
viscosity = 0.00214  # ---------------------------抛光液粘度
k = hpiv**3 / r0 * p0 / r0
xx = r0 / hpiv
aa = 6 * viscosity * wp / p0 * xx**2
dd = d / r0
ee = ww / wp
rho = 1800.0

# Create mesh and define function space
domain = Circle(Point(0, 0), 1)
mesh = generate_mesh(domain, 64)
V = FunctionSpace(mesh, 'P', 2)

# Define boundary condition
w_D = Constant(1.0)


def boundary(x, on_boundary):
    return on_boundary


bc = DirichletBC(V, w_D, boundary)

# Define load
p = Expression('1-kk*x[0]*sin(angle1)-kk*x[1]*sin(angle2)',
               degree=2,
               kk=r0 / hpiv,
               angle1=angle1,
               angle2=angle2)

u1_0 = Expression('wp*(x[0]*r0+d)', degree=2, wp=wp, d=d, r0=r0)
u1_h = Expression('ww*x[0]*r0', degree=2, ww=ww, r0=r0)
u2_0 = Expression('-wp*x[1]*r0', degree=2, wp=wp, r0=r0)
u2_h = Expression('-ww*x[1]*r0', degree=2, ww=ww, r0=r0)
dx_u1_0 = wp * r0
dx_u1_h = ww * r0
dy_u2_0 = -wp * r0
dy_u2_h = -ww * r0
f = Expression('rho*(-x[0]*r0*ww*ww)', degree=2, rho=rho, ww=ww, r0=r0)
g = Expression('rho*(-x[1]*r0*ww*ww)', degree=2, rho=rho, ww=ww, r0=r0)
px_p = -xx * sin(angle1)
py_p = -xx * sin(angle2)
ppp = Expression('p*p*p', degree=2, p=p)
fun1 = Expression('r0*r0/p0*(ww*ww*ppp+3*f*p*p/r0*px_p)',
                  degree=2,
                  ppp=ppp,
                  p=p,
                  r0=r0,
                  f=f,
                  px_p=px_p,
                  p0=p0,
                  ww=ww)
fun2 = Expression('r0*r0/p0*(ww*ww*ppp+3*g*p*p/r0*py_p)',
                  degree=2,
                  ppp=ppp,
                  p=p,
                  r0=r0,
                  g=g,
                  py_p=py_p,
                  p0=p0,
                  ww=ww)
F0 = Expression('rho*u1_0*wp', degree=2, rho=rho, u1_0=u1_0, wp=wp)
Fh = Expression('rho*u1_h*ww', degree=2, rho=rho, u1_h=u1_h, ww=ww)
G0 = Expression('-rho*u2_0*wp', degree=2, rho=rho, u2_0=u2_0, wp=wp)
Gh = Expression('-rho*u2_h*ww', degree=2, rho=rho, u2_h=u2_h, ww=ww)
F = Expression('Fh - 2*F0', degree=2, Fh=Fh, F0=F0)
G = Expression('Gh - 2*G0', degree=2, Gh=Gh, G0=G0)
d_Fx = rho * wp**2 * (ee**2 - 2) * r0
d_Gy = d_Fx
duiliu = Expression('r0/p0*(ppp*d_Fx+3*p*p*F*px_p+ppp*d_Gy+3*p*p*G*py_p)',
                    degree=2,
                    ppp=ppp,
                    d_Fx=d_Fx,
                    p=p,
                    F=F,
                    px_p=px_p,
                    d_Gy=d_Gy,
                    G=G,
                    py_p=py_p,
                    r0=r0,
                    p0=p0)

f = Expression(
    '6*viscosity*((wp*(x[0]*r0+d)+ww*r0*x[0])*sin(angle1)'
    '-((ww+wp)*r0*x[1])*sin(angle2))/k-fun1-fun2',
    degree=2,
    viscosity=viscosity,
    wp=wp,
    ww=ww,
    r0=r0,
    angle1=angle1,
    angle2=angle2,
    k=k,
    d=d,
    fun1=fun1,
    fun2=fun2)

# Define variational problem
w = TrialFunction(V)
v = TestFunction(V)
# a = dot(grad(w), grad(v))*dx
a = p * inner(grad(w), grad(v)) * dx
L = f * v * dx

# Compute solution
w = Function(V)
solve(a == L, w, bc)
plot(w, title='Deflection', cmap=cm.rainbow)
print(np.max(np.array(w.vector())))
print(np.min(np.array(w.vector())))
# Save solution to file in VTK format
# vtkfile_w = File('cmpdata/duiliu.pvd')
# vtkfile_w << w

import matplotlib.pyplot as plt

plt.show()
