import os

import matplotlib

if os.name == "posix" and 'DISPLAY' not in os.environ:
    print("MATPLOTLIB: No Display found, using non-interactive svg backend", sys.stderr)
    matplotlib.use('svg')
    import matplotlib.pyplot as plt
else:
    # matplotlib.use('svg')
    import matplotlib.pyplot as plt

    plt.ion()

__author__ = 'fzhang'

TEST_MTPY_ROOT = os.path.normpath(
    os.path.abspath(
        os.path.dirname(
            os.path.dirname(__file__)
        )
    )
)  # assume tests is on the root level of mtpy

TEST_DIR = os.path.normpath(os.path.abspath(os.path.dirname(__file__)))
TEST_TEMP_DIR = os.path.normpath(os.path.join(TEST_DIR, "temp"))

if not os.path.isdir(TEST_TEMP_DIR):
    os.mkdir(TEST_TEMP_DIR)


def _plt_wait(seconds):
    if plt.isinteractive() and plt.get_fignums():
        plt.pause(seconds)


def _plt_close():
    if plt.get_fignums():
        plt.close("all")
