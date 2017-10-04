from __future__ import print_function

import os
import matplotlib

import sys

if os.name == "posix" and 'DISPLAY' not in os.environ:
    print("MATPLOTLIB: No Display found, using non-interactive agg backend", sys.stderr)
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
else:
    # matplotlib.use('svg')
    import matplotlib.pyplot as plt
    plt.ion()

