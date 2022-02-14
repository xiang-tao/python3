import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')


ary = np.genfromtxt("/home/xt/github/python3/cmpdata/beta_www64.txt", dtype=None)
beta = ary[0]
wx = ary[1]
wy = ary[2]
wf = ary[3]

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.grid()
lns1 = ax1.plot(beta, wf, c="b", ls="-", marker='o', label='$\overline{W}_f$')
ax1.set_xlabel("$\\beta$", fontsize='20')
ax1.set_ylabel("$\overline{W}_f$", color="b", fontsize='20')
ax1.tick_params("y", colors="b")

ax2=ax1.twinx()
lns2 = ax2.plot(beta, wx, c="r", ls="-", marker='s', label='$\overline{W}_x$')
lns3= ax2.plot(beta, wy, c="k", ls="-", marker='p', label='$\overline{W}_y$')
ax2.set_ylabel("$\overline{W}_x,\overline{W}_y$", color="r", fontsize='20')
ax2.tick_params("y", colors="r")

lns = lns1+lns2 +lns3
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc=9)

plt.show()