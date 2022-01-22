import Mesh as mesh
import LagrangeSpace as space
import Assemble as am
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
import scipy.sparse.linalg
from matplotlib import cm


point1 = mesh.Point(-1, -1)
point2 = mesh.Point(1, 1)
N1 = 64
N2 = 64
m = mesh.RectangleMesh(point1, point2, N1, N2)
# print(m.generate_boundaryedges())
V = space.LagrangeSpace(m, 1)
u = space.TrialFunction(V)
v = space.TestFunction(V)
grad_u = space.NablaTrialFunction(u)
grad_v = space.NablaTestFunction(v)


def f(x, y):
    return -y*(1-y)*(1-x-x**2/2)*np.exp(x+y)


assemble = am.Am(V)
A = assemble.assemble_matrix(grad_u, grad_v, 1, 0)
b = assemble.assemble_vector(f, v, 1)
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
# print(A)
A = sp.csr_matrix(A)
val = sp.linalg.spsolve(A, b)
print(max(val))
B = np.zeros((N2+1, N1+1))
k = 0
for i in range(N2+1):
    for j in range(N1+1):
        B[j][i] = val[k]
        k += 1
X = np.linspace(-1, 1, N1+1)
Y = np.linspace(-1, 1, N2+1)
X, Y = np.meshgrid(X, Y)
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(X, Y, B, cmap=cm.rainbow)
plt.show()



