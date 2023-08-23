import numpy as np
import matplotlib.pyplot as plt

arr = np.genfromtxt("/home/xt/github/python3/fenics/kearrdata.txt", dtype=None)
arr1 = arr[0]
keshiarr = arr[1]
odearr = np.genfromtxt("/home/xt/github/python3/fenics/odearrdata.txt",
                       dtype=None)
plt.plot(arr1, keshiarr)
plt.plot(arr1, odearr)
plt.show()
