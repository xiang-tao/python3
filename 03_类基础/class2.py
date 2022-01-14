import numpy as np

a = np.array([1,2,3,4,5]).reshape(-1,1)
print(a[[1,2,3],:][1])