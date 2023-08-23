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
angle1 = 0.015 * pi / 180  # ---------------------------转角
angle2 = 0.015 * pi / 180  # ---------------------------倾角
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
# 磁场参数
a = r0
b = r0
ahat = 1
bhat = 1
A = np.sqrt(ahat**2 + bhat**2)
I = 40
H0 = I / (2 * pi) * 1 / (np.sqrt(a**2 + b**2))
mu_0_kapai = 20
# Create mesh and define function space
domain = Circle(Point(0, 0), 1)
mesh = generate_mesh(domain, 64)
V = FunctionSpace(mesh, 'P', 2)

# Define boundary condition
w_D = Constant(1.0)
# w_D = Constant(0.0)
p = Expression('1-kk*x[0]*sin(angle1)-kk*x[1]*sin(angle2)',
               degree=2,
               kk=r0 / hpiv,
               angle1=angle1,
               angle2=angle2)
B = Expression('sqrt((x[0]-ahat)*(x[0]-ahat)+(x[1]-bhat)*(x[1]-bhat))',
               degree=2,
               ahat=ahat,
               bhat=bhat)

Hhat = Expression('A*1/B', degree=2, A=A, B=B)
Hhat_dx = Expression('-(x[0]-ahat)/(B*B*B)', degree=2, ahat=ahat, B=B)
Hhat_dy = Expression('-(x[1]-bhat)/(B*B*B)', degree=2, bhat=bhat, B=B)
Hhat_dxx = Expression(
    '(2*(x[0]-ahat)*(x[0]-ahat)-(x[1]-bhat)*(x[1]-bhat))/pow(B,5)',
    degree=2,
    B=B,
    ahat=ahat,
    bhat=bhat)
Hhat_dyy = Expression(
    '(2*(x[1]-bhat)*(x[1]-bhat)-(x[0]-ahat)*(x[0]-ahat))/pow(B,5)',
    degree=2,
    B=B,
    ahat=ahat,
    bhat=bhat)

cichang_x = Expression(
    '3*H0*H0/p0*p*p*Hhat*(-xx*sin(angle1))*Hhat_dx+H0*H0/p0*p*p*p*(Hhat_dx*Hhat_dx+Hhat*Hhat_dxx)',
    degree=2,
    p0=p0,
    H0=H0,
    xx=xx,
    angle1=angle1,
    Hhat=Hhat,
    Hhat_dx=Hhat_dx,
    Hhat_dxx=Hhat_dxx,
    p=p)
cichang_y = Expression(
    '3*H0*H0/p0*p*p*Hhat*(-xx*sin(angle2))*Hhat_dy+H0*H0/p0*p*p*p*(Hhat_dy*Hhat_dy+Hhat*Hhat_dyy)',
    degree=2,
    p0=p0,
    H0=H0,
    xx=xx,
    angle2=angle2,
    Hhat=Hhat,
    Hhat_dy=Hhat_dy,
    Hhat_dyy=Hhat_dyy,
    p=p)


def boundary(x, on_boundary):
    return on_boundary


bc = DirichletBC(V, w_D, boundary)

# Define load
# +12*viscosity*(ww*r0*x[1]*sin(angle1)-ww*r0*x[0]*sin(angle2))/k'
f = Expression(
    '6*viscosity*((wp*(x[0]*r0+d)-ww*r0*x[0])*sin(angle2)'
    '+((ww-wp)*r0*x[1])*sin(angle1))/k-mu_0_kapai*(cichang_y+cichang_x)',
    degree=2,
    mu_0_kapai=mu_0_kapai,
    cichang_x=cichang_x,
    cichang_y=cichang_y,
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
vtkfile_w = File('cmpdata/new_cichang.pvd')
vtkfile_w << w

import matplotlib.pyplot as plt

plt.show()
