#!/usr/bin/env python
import matplotlib
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, polar=True)

# make all major grid lines lighter, only x grid lines set in 2.1.0
ax.grid(alpha=0.2)   

# hide y tick labels, no effect in 2.1.0
plt.setp(ax.yaxis.get_ticklabels(), visible=False) 
fig.show()

fig.savefig("polar.png")

