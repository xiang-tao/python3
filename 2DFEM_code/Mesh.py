import numpy as np


class Point(object):
    def __init__(self, rn, cn):
        self.rn = rn
        self.cn = cn


class RectangleMesh(object):
    def __init__(self, point1, point2, N1, N2):
        self.point1 = point1
        self.point2 = point2
        self.N1 = N1
        self.N2 = N2
        self.N = N1*2*N2

    def generate_p(self):
        P = np.zeros(((self.N1 + 1)*(self.N2 + 1), 2))
        h1 = (self.point2.rn - self.point1.rn)/self.N1
        h2 = (self.point2.cn - self.point1.cn) / self.N2
        for rn in range(self.N2+1):
            for cn in range(self.N1+1):
                j = rn + cn*(self.N2+1)
                P[j][0] = self.point1.rn+cn*h1
                P[j][1] = self.point1.cn+rn*h2
        return P

    def generate_t(self):
        T = np.zeros((2*self.N1*self.N2, 3), dtype=np.int64)
        for re in range(self.N2):
            for ce in range(self.N1):
                n1 = ce*2*self.N2+2*re
                node1 = re + ce * (self.N2 + 1)
                node2 = re + (ce+1) * (self.N2 + 1)
                node3 = re + ce * (self.N2 + 1)+1
                T[n1][0] = node1
                T[n1][1] = node2
                T[n1][2] = node3
                n2 = ce*2*self.N2+2*re+1
                node1 = re + ce * (self.N2 + 1) + 1
                node2 = re + (ce + 1) * (self.N2 + 1)
                node3 = re + (ce + 1) * (self.N2 + 1)+1
                T[n2][0] = node1
                T[n2][1] = node2
                T[n2][2] = node3
        return T

    def generate_boundaryedges(self, type=-1):
        T = self.generate_t()
        boundaryedges = np.zeros((2*(self.N1+self.N2), 4))
        for ce in range(self.N1):
            boundaryedges[ce][0] = type
            nk = ce*2*self.N2
            boundaryedges[ce][1] = nk
            boundaryedges[ce][2] = T[nk][0]
            boundaryedges[ce][3] = T[nk][1]
        for re in range(self.N2):
            boundaryedges[self.N1+re][0] = type
            nk = 2*self.N2*(self.N1-1)+2*re+1
            boundaryedges[self.N1+re][1] = nk
            boundaryedges[self.N1+re][2] = T[nk][1]
            boundaryedges[self.N1+re][3] = T[nk][2]
        for ce in range(self.N1):
            boundaryedges[2*self.N1+self.N2-ce-1][0] = type
            nk = (ce+1)*2*self.N2-1
            boundaryedges[2*self.N1+self.N2-ce-1][1] = nk
            boundaryedges[2*self.N1+self.N2-ce-1][2] = T[nk][2]
            boundaryedges[2*self.N1+self.N2-ce-1][3] = T[nk][0]
        for re in range(self.N2):
            boundaryedges[2*(self.N1+self.N2)-re-1][0] = type
            nk = 2*re
            boundaryedges[2*(self.N1+self.N2)-re-1][1] = nk
            boundaryedges[2*(self.N1+self.N2)-re-1][2] = T[nk][2]
            boundaryedges[2*(self.N1+self.N2)-re-1][3] = T[nk][0]
        return boundaryedges
