import numpy as np
import frogpy.algebra.matrix as matrix
a = matrix.tri_diagonal(4, 2, 2, 1)
print(a)


def test(a, b):
    """
    :param a:integer
    :param b:integer
    """
    c = a+b
    print(c)


a = 3
b = 4
test(a, b)

print(np.pi)
print("hello")
