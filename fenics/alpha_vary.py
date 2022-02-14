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
hpiv = 8.0*1e-5  # ----------------------------------------晶片中心高度
viscosity = 0.00214  # ---------------------------抛光液粘度
k = hpiv**3/r0*p0/r0
xx = r0 / hpiv
aa = 6 * viscosity * wp / p0 * xx ** 2
dd = d / r0
ee = ww / wp
rho = 1800.0
num = 7
alpha_arry = np.linspace(-0.03, 0.03, num)
wf = np.zeros(num)
wx = np.zeros(num)
wy = np.zeros(num)


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


def cell_values_avg_px(coorcell, u):
    n = len(coorcell)
    avg = 0.0
    for i in range(n):
        avg += u(coorcell[i])*coorcell[i][1]
    return avg/n


def cell_values_avg_py(coorcell, u):
    n = len(coorcell)
    avg = 0.0
    for i in range(n):
        avg += u(coorcell[i])*coorcell[i][0]
    return -avg/n


for j in range(num):
    angle1 = alpha_arry[j] * pi / 180
    k = hpiv ** 3 / r0 * p0 / r0
    xx = r0 / hpiv
    aa = 6 * viscosity * wp / p0 * xx ** 2
    dd = d / r0
    ee = ww / wp
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
                   degree=2, kk=r0 / hpiv, angle1=angle1, angle2=angle2)

    u1_0 = Expression('wp*(x[0]*r0+d)', degree=2, wp=wp, d=d, r0=r0)
    u1_h = Expression('ww*x[0]*r0', degree=2, ww=ww, r0=r0)
    u2_0 = Expression('-wp*x[1]*r0', degree=2, wp=wp, r0=r0)
    u2_h = Expression('-ww*x[1]*r0', degree=2, ww=ww, r0=r0)
    dx_u1_0 = wp * r0
    dx_u1_h = ww * r0
    dy_u2_0 = -wp * r0
    dy_u2_h = -ww * r0
    F0 = Expression('rho*u1_0*wp', degree=2, rho=rho, u1_0=u1_0, wp=wp)
    Fh = Expression('rho*u1_h*ww', degree=2, rho=rho, u1_h=u1_h, ww=ww)
    G0 = Expression('-rho*u2_0*wp', degree=2, rho=rho, u2_0=u2_0, wp=wp)
    Gh = Expression('-rho*u2_h*ww', degree=2, rho=rho, u2_h=u2_h, ww=ww)
    F = Expression('Fh - 2*F0', degree=2, Fh=Fh, F0=F0)
    G = Expression('Gh - 2*G0', degree=2, Gh=Gh, G0=G0)
    ppp = Expression('p*p*p', degree=2, p=p)
    px_p = -xx * sin(angle1)
    py_p = -xx * sin(angle2)
    d_Fx = rho * wp ** 2 * (ee ** 2 - 2) * r0
    d_Gy = d_Fx
    duiliu = Expression('r0/p0*(ppp*d_Fx+3*p*p*F*px_p+ppp*d_Gy+3*p*p*G*py_p)',
                        degree=2, ppp=ppp, d_Fx=d_Fx, p=p, F=F, px_p=px_p,
                        d_Gy=d_Gy, G=G, py_p=py_p, r0=r0, p0=p0)

    f = Expression('6*viscosity*((wp*(x[0]*r0+d)+ww*r0*x[0])*sin(angle1)'
                   '-((ww+wp)*r0*x[1])*sin(angle2))/k-duiliu',
                   degree=2, viscosity=viscosity, wp=wp, ww=ww, r0=r0,
                   angle1=angle1, angle2=angle2, k=k, d=d, duiliu=duiliu)

    # Define variational problem
    w = TrialFunction(V)
    v = TestFunction(V)
    # a = dot(grad(w), grad(v))*dx
    a = p * inner(grad(w), grad(v)) * dx
    L = f * v * dx

    # Compute solution
    w = Function(V)
    solve(a == L, w, bc)
    # plot(w, title='Deflection', cmap=cm.rainbow)

    load = 0.0
    int_px = 0.0
    int_py = 0.0
    element = V.element()
    for cell in cells(mesh):
        coorcell = element.tabulate_dof_coordinates(cell)
        measure = cells_measure(coorcell)
        values_avg = cell_values_avg(coorcell, w)
        load += values_avg*measure
        int_px += cell_values_avg_px(coorcell, w)*measure
        int_py += cell_values_avg_py(coorcell, w)*measure

    wf[j] = load/pi
    wx[j] = int_px/pi
    wy[j] = int_py/pi

cmpwww = np.array([alpha_arry, wx, wy, wf])

np.savetxt('/home/xt/github/python3/fenics/cmpdata/alpha_vary_data.txt', np.c_[cmpwww],
           fmt='%.16f', delimiter='\t')