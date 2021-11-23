import numpy as np
from fealpy.decorator import cartesian


class PdeData:

    def __init__(self):
        pass

    def domain(self):
        return np.array([0, 1, 0, 1])

    @cartesian
    def solution(self, p):
        """ The exact solution
        Parameters
        ---------
        p :


        Examples
        -------
        p = np.array([0, 1], dtype=np.float64)
        p = np.array([[0, 1], [0.5, 0.5]], dtype=np.float64)
        """
        x = p[..., 0]
        y = p[..., 1]
        val = x**2+y**2-1
        return val  # val.shape == x.shape

    @cartesian
    def source(self, p):
        x = p[..., 0]
        y = p[..., 1]
        val = -6*(x+y)

        return val

    @cartesian
    def dirichlet(self, p):
        """
        边界处理
        """
        return self.solution(p)
        # return 0

    @cartesian
    def h_function(self, p):
        x = p[..., 0]
        y = p[..., 1]
        val = x+y
        return val  # val.shape == x.shape
