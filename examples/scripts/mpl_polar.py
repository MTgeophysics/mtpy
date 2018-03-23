#!/usr/bin/env python
"""
This script demonstrates how to make a polar plot

Author: unknown

"""


import matplotlib
import matplotlib.pyplot as plt

print(matplotlib.__version__)  # print version
print(matplotlib.get_backend())  # print backend

fig = plt.figure()
ax = fig.add_subplot(111, polar=True)

# make all major grid lines lighter, only x grid lines set in 2.1.0
ax.grid(alpha=0.2)   

# hide y tick labels, no effect in 2.1.0
plt.setp(ax.yaxis.get_ticklabels(), visible=False) 
fig.show()

fig.savefig("polar.png")

# prints alpha of one major gridline of each axis
print(ax.xaxis.majorTicks[0].gridline._alpha)
print(ax.yaxis.majorTicks[0].gridline._alpha)
