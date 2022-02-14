import fem2 as fem
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

point1 = fem.Point(0, 0)
point2 = fem.Point(1, 1)
N1 = 16
N2 = 16
m = fem.RectangleMesh(point1, point2, N1, N2)
V = fem.FunctionSpace(m, degree=1)
u = fem.TrialFunction(V)
v = fem.TestFunction(V)
grad_u = fem.nabla_u(u)
grad_v = fem.nabla_v(v)


def f(x, y):
    return 2*np.pi**2*np.sin(np.pi*x)*np.sin(np.pi*y)


def solution(x, y):
    return np.sin(np.pi*x)*np.sin(np.pi*y)


A = fem.assemble_matrix(V, grad_u, grad_v)
b = fem.assemble_vector(V, f, v)
boundarynodes = V.generate_boundarynodes()
Pb = V.generate_pb()
for i in range(boundarynodes.shape[0]):
    if boundarynodes[i][0] == -1:
        A[boundarynodes[i][1], :] = 0
        A[boundarynodes[i][1]][boundarynodes[i][1]] = 1
        b[boundarynodes[i][1]] = 0

uh = fem.Function(V)
fem.solve(A, uh, b)
L2error = fem.L2Error(solution, uh)
print(L2error)
val = uh.Vector()
print(max(val))
B = np.zeros((N2+1, N1+1))
k = 0
for i in range(N2+1):
    for j in range(N1+1):
        B[j][i] = val[k]
        k += 1
X = np.linspace(0, 1, N1+1)
Y = np.linspace(0, 1, N2+1)
X, Y = np.meshgrid(X, Y)
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(X, Y, B, cmap=cm.rainbow)
plt.show()