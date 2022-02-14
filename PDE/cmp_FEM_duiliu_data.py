import numpy as np
from fealpy.decorator import cartesian

sin = np.sin
cos = np.cos
pi = np.pi
# 请注意r/s单位不是国际单位，应换算成弧度每秒，而一转的弧度是2*pi，故需要乘2*pi
ww = 5.0 * pi / 3.0  # 50.0  # ----------------------------------------晶片的角速度
wp = 10.0 * pi / 3.0  # 100.0  # ----------------------------------------抛光垫的角速度
angle1 = 0.02 * pi / 180  # ---------------------------转角
angle2 = 0.018 * pi / 180  # ---------------------------倾角
d = 0.15  # ---------------------------晶片和抛光垫的旋转中心距
r0 = 5.0*1e-2
p0 = 101000.0
hpiv = 8.0*1e-5  # ----------------------------------------晶片中心高度
viscosity = 0.00214  # ---------------------------抛光液粘度
k = hpiv**3/r0*p0/r0
xx = r0 / hpiv
aa = 6 * viscosity * wp / p0 * xx ** 2
dd = d / r0
ee = ww / wp
rho = 1800.0


def h(self, p):
    x = p[..., 0]
    y = p[..., 1]
    return 1-xx*x*sin(angle1)-xx*y*sin(angle2)


def dh_x(self, p):
    x = p[..., 0]
    y = p[..., 1]
    return -xx*sin(angle1)


def dh_y(self, p):
    x = p[..., 0]
    y = p[..., 1]
    return -xx*sin(angle2)


def B(self, p):
    x = p[..., 0]
    return rho*ee**2*wp**2*r0*x-2*rho*wp**2*(r0*x+d)


dB_x = rho*wp**2*(ee**2-2)*r0


def C(self, p):
    y = p[..., 1]
    return rho*ee**2*wp**2*r0*y-2*rho*wp**2*r0*y


dC_y = dB_x


def duiliu(self, p):
    x = p[..., 0]
    y = p[..., 1]
    return r0/p0*(h(self, p)**3*dB_x+B(self, p)*3*h(self, p)**2
                  *dh_x(self, p)+h(self, p)**3*dC_y+C(self, p)*3*h(self, p)**2*dh_y(self, p))


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

        val = 6*viscosity*((wp*(x*r0+d)+ww*r0*x)*sin(angle1)-((ww+wp)*r0*y)*sin(angle2))/k-duiliu(self, p)
        return val

    @cartesian
    def dirichlet(self, p):
        """
        边界处理
        """
        pp = 1.0  # ---------------------------无量纲后的大气压边界初始值
        return pp

    @cartesian
    def h_function(self, p):
        """
        定义厚度模型h
        """
        x = p[..., 0]
        y = p[..., 1]
        kk = r0/hpiv
        val = 1-kk*x*sin(angle1)-kk*y*sin(angle2)
        return val  # val.shape == x.shape
