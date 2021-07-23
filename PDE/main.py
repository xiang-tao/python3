import numpy as np
from fealpy.mesh import TriangleMesh
import matplotlib.pyplot as plt

B = np.ones((3, 3))
L = [2,2,3]
B[:,0] = L
v0 = [1,2,3]
print(L[1])
