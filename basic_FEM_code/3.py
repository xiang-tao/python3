import numpy as np
import Mesh, Space, Assemble, BoundaryProcess

left, right, N = 0, 1, 5
mesh = Mesh.IntervalMesh(left, right, N)
print(mesh.generate_p())