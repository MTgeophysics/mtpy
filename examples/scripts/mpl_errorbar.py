"""
This script demonstrates how to make a plot with errorbars that go below zero 
on a semilog plot

Author: unknown

"""

import matplotlib
import matplotlib.pyplot as plt

print(matplotlib.__version__)
print(matplotlib.get_backend())

x = range(100)
y = [50] * 100
yerr = range(100)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.errorbar(x, y, yerr=yerr)

ax.set_yscale("log", nonposy="clip")
# ax.set_xscale('log')
# ax.set_ylim((10, 100))

# fig.savefig("errorbar_mpl{}.png".format(matplotlib.__version__.replace(".", "")))

plt.show()
