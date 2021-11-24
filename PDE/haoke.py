import xt_pdenew_matrix as matrix
import numpy as np
a = 1.58
b = 1.2
nx = 100
ny = 100
A = matrix.tri_diagonal(nx, 1+2*a, -a, -a)
arr1 = np.ones(nx)
arr2 = np.ones(nx)*-1
B = np.zeros((nx,nx))
matrix.fill_diag(B,arr1,0)
matrix.fill_diag(B,arr2,-1)
C = B *b
D = B*-b
ex = np.zeros(100)
x = np.linspace(0,1,100)
y = np.linspace(0,1,100)
z = np.linspace(0,1,100)

def f(x, y, z):
    return x+y+z

for k in range(100):
    ex[k] = f(1,1,z[k])
C = B@ex
print(C)
print(C.shape)