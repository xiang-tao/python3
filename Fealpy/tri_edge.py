import numpy as np
import matplotlib.pyplot as plt
from fealpy.mesh import TriangleMesh

node = np.array((
    [0.0, 0.0],
    [1.0, 0.0],
    [1.0, 1.0],
    [0.0, 1.0]
), dtype=np.float64)