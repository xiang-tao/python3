from fenics import *
import numpy as np
import matplotlib.pyplot as plt


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


mesh = UnitSquareMesh(1, 1)
V = FunctionSpace(mesh, 'P', 2)
u = interpolate(Expression('x[0]+x[1]', degree=1), V)
# plot(u)
# nodal_values = np.array(u.vector())
# vertex_values = np.array(u.compute_vertex_values(mesh))
# print(nodal_values)
# print(vertex_values)
element = V.element()
for cell in cells(mesh):
    # print(element.tabulate_dof_coordinates(cell))
    coorcell = element.tabulate_dof_coordinates(cell)
    measure = cells_measure(coorcell)
    print(len(coorcell))
    print(measure)

# coordinates = mesh.coordinates()
# print(coordinates)

# plt.show()