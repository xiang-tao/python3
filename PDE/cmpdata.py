import numpy as np
from fealpy.decorator import cartesian


class CmpData:

    def __init__(self):
        pass

    def domain(self):
        return np.array([0, 1, 0, 1])

    @cartesian
    def source(self, p):
        """
        定义求解右端载荷向量函数
        定义PDE模型，给定无量纲后的数据初始值
        """
        x = p[..., 0]
        y = p[..., 1]
        pi = np.pi
        ww = 5.0 * pi / 3.0  # 50.0  # ----------------------------------------晶片的角速度
        wp = 10.0 * pi / 3.0  # 100.0  # ----------------------------------------抛光垫的角速度
        angle1 = 0.02 * pi / 180  # ---------------------------转角
        angle2 = 0.018 * pi / 180  # ---------------------------倾角
        d = 150000.0  # ---------------------------晶片和抛光垫的旋转中心距
        r0 = 50000.0
        p0 = 101000
        hpiv = 80.0  # ----------------------------------------晶片中心高度
        viscosity = 0.00214  # ---------------------------抛光液粘度
        k = hpiv**3/r0*p0/r0
        val = 6*viscosity*((wp*(x*r0+d)+ww*r0*x)*np.sin(angle1)-((ww+wp)*r0*y)*np.sin(angle2))/k
        return val

    @cartesian
    def dirichlet(self, p):
        """
        边界处理
        """
        p0 = 1.0  # ---------------------------无量纲后的大气压边界初始值
        return p0

    @cartesian
    def h_function(self, p):
        """
        定义厚度模型h
        """
        x = p[..., 0]
        y = p[..., 1]
        hpiv = 80.0
        pi = np.pi
        angle1 = 0.02 * pi / 180  # ---------------------------转角
        angle2 = 0.018 * pi / 180  # ---------------------------倾角
        r0 = 50000.0
        kk = r0/hpiv
        # val = hpiv-x*np.sin(angle1)-y*np.sin(angle2)
        val = 1-kk*x*np.sin(angle1)-kk*y*np.sin(angle2)
        return val  # val.shape == x.shape
