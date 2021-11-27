import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')


ary = np.genfromtxt("/home/xt/github/python3/cmpdata/wfwxwy.txt", dtype=None)
hpiv = ary[0]
wx = ary[1]
wy = ary[2]
wf = ary[3]

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.grid()
lns1 = ax1.plot(hpiv, wf, c="b", ls="-", marker='o', label='wf')
ax1.set_xlabel("x")
ax1.set_ylabel("e", color="b")
ax1.tick_params("y", colors="b")

ax2=ax1.twinx()
lns2 = ax2.plot(hpiv, wx, c="r", ls="-", marker='s', label='wx')
lns3= ax2.plot(hpiv, wy, c="k", ls="-", marker='p', label='wy')
# ax2.set_ylabel("yuxuan", color="r")
# ax2.tick_params("y", colors="r")

lns = lns1+lns2+lns3
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc=0)

plt.show()