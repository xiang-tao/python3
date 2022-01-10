import numpy as np
import scipy.sparse as sp
from scipy import integrate
import scipy.sparse.linalg
import matplotlib.pyplot as plt
from matplotlib import cm

def solve():
    P, T = generate_PbTb(?,101)
    if basis_type_trial==101:
        Pb_trial =
        Tb_trial =
    elif basis_type_trial==102:
        pass
    else:
        pass

    if basis_type_test==101:
        Pb_test =
        Tb_test =
    elif basis_type_test==102:
        pass
    else:
        pass


boundarynodes = generate_boundarynodes(?)
A = assemble_matrix(?)
b = assemble_vector(?)

A, b = DirichletBC(?)
