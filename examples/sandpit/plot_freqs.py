"""
Description:
    what does this script module do? How to do it.

Author: fei.zhang@ga.gov.au

Date:
"""

__author__ = 'u25656'

def plot_freqlist():

    # import matplotlib
    # matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt

    freqs = [ 1.000000e+00,  5.623413e-01, 3.162278e-01 , 1.778279e-01 , 1.000000e-01 ,5.623413e-02,
     3.162278e-02 , 1.778279e-02 , 1.000000e-02 , 5.623413e-03 , 3.162278e-03 , 1.778279e-03,
     1.000000e-03,  5.623413e-04 , 3.162278e-04 , 1.778279e-04 , 1.000000e-04 ]

    periods= [1/f for f in freqs]
    plt.plot(freqs, '^', markersize='10')
    plt.show()

    plt.plot(periods, 'o', markersize='10')
    plt.show()

def try_scipy_interp():

    import numpy as np
    from scipy.interpolate import interp1d
    import matplotlib.pyplot as plt

    x = [-1.01, 5.66, 5.69, 13.77, 20.89]
    y = [0.28773, 1.036889, 1.043178, 1.595322, 1.543763]

    new_x = [0, 2, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20]
    new_y_scipy = [0.401171, 0.625806, 0.850442, 1.062384, 1.186291, 1.248244, 1.310198, 1.372152, 1.434105, 1.496059,
                 1.545429, 1.55267, 1.559911, 1.567153, 1.574394, 1.588877, ]
    new_y_matlab = [0.401171, 0.625806, 0.850442, 1.064362, 1.201031, 1.269366, 1.3377, 1.406035, 1.47437, 1.542704,
                  1.593656, 1.586415, 1.579174, 1.571932, 1.564691, 1.550208]

    askewchan = interp1d(x, y)(new_x)

    print askewchan

    # 'linear' has no effect since it's the default, but I'll plot it too:
    set_interp = interp1d(x, y, kind='linear')
    new_y = set_interp(new_x)

    print new_y

    # plt.plot(x, y, 'o', new_x, new_y_scipy, '--', new_x, new_y_matlab, ':', new_x, askewchan, '^', new_x, new_y, '+')
    # plt.legend(('Original', 'OP_scipy', 'OP_matlab', 'askewchan_scipy', 'OP style scipy'), loc='lower right')
    plt.plot(x, y, 'o', new_x, askewchan, '^', new_x, new_y, '+')
    plt.legend(('Original', 'OP_scipy', 'OP_matlab', 'askewchan_scipy', 'OP style scipy'), loc='lower right')

    np.allclose(new_y_matlab, interp1d(x, y)(new_x))
    # True
    plt.show()

def try_interp1D():
    # http://stackoverflow.com/questions/27698604/what-do-the-different-values-of-the-kind-argument-mean-in-scipy-interpolate-inte

    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.interpolate as interpolate

    np.random.seed(6)
    kinds = ('nearest', 'zero', 'linear', 'slinear', 'quadratic', 'cubic')

    N = 10
    x = np.linspace(0, 1, N)
    y = np.random.randint(10, size=(N,))

    new_x = np.linspace(0, 1, 28)
    fig, axs = plt.subplots(nrows=len(kinds) + 1, sharex=True)
    axs[0].plot(x, y, 'bo-')
    axs[0].set_title('raw')
    for ax, kind in zip(axs[1:], kinds):
        new_y = interpolate.interp1d(x, y, kind=kind)(new_x)
        ax.plot(new_x, new_y, 'ro-')
        ax.set_title(kind)

    plt.show()

####################################################
if __name__ == "__main__":
    #try_scipy_interp()

    try_interp1D()
