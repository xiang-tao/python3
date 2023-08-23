"""
FEniCS tutorial demo program: Deflection of a membrane.

  -Laplace(w) = p  in the unit circle
            w = 0  on the boundary

The load p is a Gaussian function centered at (0, 0.6).
"""

# from __future__ import print_function
from matplotlib import cm
# from fenics import *
from dolfin import *
from mshr import *
import numpy as np

sin = np.sin
cos = np.cos
pi = np.pi
# 请注意r/s单位不是国际单位，应换算成弧度每秒，而一转的弧度是2*pi，故需要乘2*pi
ww = 5.0 * pi / 3.0  # 50.0  # ----------------------------------------晶片的角速度
wp = 10.0 * pi / 3.0  # 100.0  # ----------------------------------------抛光垫的角速度
angle1 = 0.03 * pi / 180  # ---------------------------转角
angle2 = 0.025 * pi / 180  # ---------------------------倾角
d = 0.15  # ---------------------------晶片和抛光垫的旋转中心距
r0 = 5.0 * 1e-2
p0 = 101000.0
hpiv = 5.0 * 1e-5  # ----------------------------------------晶片中心高度
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
M = Expression('4*rho*ww*wp*(r0*x[0]+d)-2*rho*ww*ww*r0*x[0]',
               degree=2,
               r0=r0,
               rho=rho,
               ww=ww,
               wp=wp,
               d=d)
N = Expression('4*rho*ww*wp*r0*x[1]-2*rho*ww*ww*r0*x[1]',
               degree=2,
               r0=r0,
               rho=rho,
               ww=ww,
               wp=wp)
dM = 2 * rho * ww * (2 * wp - ww)
dN = dM
keshili1 = Expression('(p*p*p*dM-3*M*p*p/hpiv*sin(angle1))*r0*r0/p0',
                      degree=2,
                      p=p,
                      M=M,
                      dM=dM,
                      angle1=angle1,
                      r0=r0,
                      p0=p0,
                      hpiv=hpiv)
keshili2 = Expression('(p*p*p*dN-3*N*p*p/hpiv*sin(angle2))*r0*r0/p0',
                      degree=2,
                      p=p,
                      N=N,
                      dN=dN,
                      angle2=angle2,
                      r0=r0,
                      p0=p0,
                      hpiv=hpiv)
f = Expression(
    '6*viscosity*((wp*(x[0]*r0+d)-ww*r0*x[0])*sin(angle2)'
    '+((ww-wp)*r0*x[1])*sin(angle1))/k-keshili2-keshili1',
    keshili1=keshili1,
    keshili2=keshili2,
    degree=2,
    viscosity=viscosity,
    wp=wp,
    ww=ww,
    r0=r0,
    angle1=angle1,
    angle2=angle2,
    k=k,
    d=d)

# Define variational problem
w = TrialFunction(V)
v = TestFunction(V)
# a = dot(grad(w), grad(v))*dx
a = p * inner(grad(w), grad(v)) * dx  # type:ignore
L = f * v * dx  # type:ignore
# Compute solution
w = Function(V)
solve(a == L, w, bc)
plot(w, title='Deflection', cmap=cm.rainbow)  # type:ignore
print(np.max(np.array(w.vector())))
print(np.min(np.array(w.vector())))

# Save solution to file in VTK format
# vtkfile_w = File('cmpdata/new_keshili.pvd')
# vtkfile_w << w

import matplotlib.pyplot as plt

plt.show()
