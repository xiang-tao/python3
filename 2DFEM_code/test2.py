import fem2 as fem
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

point1 = fem.Point(-1, -1)
point2 = fem.Point(1, 1)
N1 = 32
N2 = 32
m = fem.RectangleMesh(point1, point2, N1, N2)
V = fem.FunctionSpace(m, degree=1)
u = fem.TrialFunction(V)
v = fem.TestFunction(V)
grad_u = fem.nabla_u(u)
grad_v = fem.nabla_v(v)


def f(x, y):
    return -y*(1-y)*(1-x-x**2/2)*np.exp(x+y)-x*(1-x/2)*(-3*y-y**2)*np.exp(x+y)


def c(x, y):
    return 1


def solution(x, y):
    return x*y*(1-x/2)*(1-y)*np.exp(x+y)


def dux(x, y):
    return np.exp(y+x)*(x-x**2/2)*(1-y)


def duy(x, y):
    return -y*np.exp(y+x)*(x-x**2/2)


A = fem.assemble_matrix(V, grad_u, grad_v)
b = fem.assemble_vector(V, f, v)
boundarynodes = V.generate_boundarynodes()
Pb = V.generate_pb()
for i in range(boundarynodes.shape[0]):
    x = Pb[boundarynodes[i][1]][0]
    y = Pb[boundarynodes[i][1]][1]
    if boundarynodes[i][0] == -1:
        A[boundarynodes[i][1], :] = 0
        A[boundarynodes[i][1]][boundarynodes[i][1]] = 1
        if x == -1:
            b[boundarynodes[i][1]] = -1.5*y*(1-y)*np.exp(-1+y)
        if x == 1:
            b[boundarynodes[i][1]] = 0.5*y*(1-y)*np.exp(1+y)
        if y == -1:
            b[boundarynodes[i][1]] = -2*x*(1-x/2)*np.exp(x-1)
        if y == 1:
            b[boundarynodes[i][1]] = 0

uh = fem.Function(V)
fem.solve(A, uh, b)
L2error = fem.L2Error(solution, uh)
H1halferror = fem.H1half(dux, duy, uh, fem.nabla_u(u))
print(H1halferror)
print(L2error)
val = uh.Vector()
print(max(val))
# B = np.zeros((N2+1, N1+1))
# k = 0
# for i in range(N2+1):
#     for j in range(N1+1):
#         B[j][i] = val[k]
#         k += 1
# X = np.linspace(-1, 1, N1+1)
# Y = np.linspace(-1, 1, N2+1)
# X, Y = np.meshgrid(X, Y)
# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# surf = ax.plot_surface(X, Y, B, cmap=cm.rainbow)
# plt.show()



