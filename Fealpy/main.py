import numpy as np
from fealpy.mesh import QuadrangleMesh
import matplotlib.pyplot as plt


def arrayprint(name, a):
    """
    Note
    ----
    打印一个名字为 name 的数组 a，每一行之前加一个行号
    """
    print("\n", name + ":")
    for (i, row) in enumerate(a):
        print(i, ": ", row)


node = np.array([(0, 0), (1, 0), (1, 1), (0, 1)], dtype=np.float)
cell = np.array([(0, 1, 2, 3)], dtype=np.int)
mesh = QuadrangleMesh(node, cell)
mesh.uniform_refine(1)

node = mesh.entity('node')
cell = mesh.entity('cell')

# 从 cell 出发，构造 edge、edge2cel、cell2edge
NC = mesh.number_of_cells()
NEC = 4  # 每个单元有4条边
localEdge = np.array([
    [0, 1],  # 局部 0 号边
    [1, 2],  # 局部 1 号边
    [2, 3],  # 局部 2 号边
    [3, 0],  # 局部 3 号边
], dtype=np.int_)  # (4, 2)
# （NC, 4)---> (NC, 4, 2) --> (4*NC, 2)
totalEdge = cell[:, localEdge].reshape(-1, 2)
stotalEdge = np.sort(totalEdge, axis=-1)

"""
这里对python的np.unique函数作说明：
1.a = np.unique(A)
对于一维数组或者列表，unique函数去除其中重复的元素，并按元素由大到小返回一个新的无元素重复的元组或者列表
2.c,s=np.unique(b,return_index=True)
return_index=True表示返回新列表元素在旧列表中的位置，并以列表形式储存在s中。
3.a, s,p = np.unique(A, return_index=True, return_inverse=True)
return_inverse=True 表示返回旧列表元素在新列表中的位置，并以列表形式储存在p中
"""
sedge, i0, j = np.unique(stotalEdge,
                         axis=0,
                         return_index=True,
                         return_inverse=True)

i1 = np.zeros_like(i0)
i1[j] = range(NEC * NC)  # 取重复元素的后者

edge = totalEdge[i0]  # 去除重复元素后将顺序复原得到边的约定数据
cell2edge = j.reshape(-1, NEC)  # (NC, 4)

NE = len(edge)
edge2cell = np.zeros((NE, 4), dtype=np.int_)
edge2cell[:, 0] = i0 // NEC
edge2cell[:, 1] = i1 // NEC
edge2cell[:, 2] = i0 % NEC
edge2cell[:, 3] = i1 % NEC

# 输出构造出的edge cell2edge edge2cell

arrayprint("edge", edge)
arrayprint("edge2cell", edge2cell)
arrayprint("cell2edge", cell2edge)

fig = plt.figure()
axes = fig.gca()
mesh.add_plot(axes)
mesh.find_node(axes, showindex=True)
mesh.find_cell(axes, showindex=True)
mesh.find_edge(axes, showindex=True)
plt.show()
