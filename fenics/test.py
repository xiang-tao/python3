"""
FEniCS tutorial demo program: Deflection of a membrane.

  -Laplace(w) = p  in the unit circle
            w = 0  on the boundary

The load p is a Gaussian function centered at (0, 0.6).
"""

# from __future__ import print_function
from matplotlib import cm
import matplotlib.pyplot as plt
from dolfin import *
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
r0 = 5.0*1e-2
p0 = 101000.0
hpiv = 5.0*1e-5  # ----------------------------------------晶片中心高度
viscosity = 0.00214  # ---------------------------抛光液粘度
k = hpiv**3/r0*p0/r0
xx = r0 / hpiv
aa = 6 * viscosity * wp / p0 * xx ** 2
dd = d / r0
ee = ww / wp
rho = 1800.0

def cells_measure(coorcell):
    coor1 = coorcell[0]
    coor2 = coorcell[1]
    coor3 = coorcell[2]
    a = np.linalg.norm(coor1-coor2)
    b = np.linalg.norm(coor1 - coor3)
    c = np.linalg.norm(coor3 - coor2)
    p = (a+b+c)/2
    ss = p*(p-a)*(p-b)*(p-c)
    s = np.sqrt(ss)
    return s


def cell_values_avg(coorcell, u):
    n = len(coorcell)
    avg = 0.0
    for i in range(n):
        avg += u(coorcell[i])
    return avg/n


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
               degree=2, kk=r0/hpiv, angle1=angle1, angle2=angle2)
f = Expression('6*viscosity*((wp*(x[0]*r0+d)+ww*r0*x[0])*sin(angle1)'
               '-((ww+wp)*r0*x[1])*sin(angle2))/k',
               degree=2, viscosity=viscosity, wp=wp, ww=ww, r0=r0,
               angle1=angle1, angle2=angle2, k=k, d=d)

# Define variational problem
w = TrialFunction(V)
v = TestFunction(V)
# a = dot(grad(w), grad(v))*dx
a = p*inner(grad(w), grad(v))*dx
L = f*v*dx

# Compute solution
w = Function(V)
solve(a == L, w, bc)
plot(w, title='Deflection', cmap=cm.rainbow)

load = 0.0
element = V.element()
for cell in cells(mesh):
    coorcell = element.tabulate_dof_coordinates(cell)
    measure = cells_measure(coorcell)
    values_avg = cell_values_avg(coorcell, w)
    load += measure*values_avg

print(load/pi)
plt.show()