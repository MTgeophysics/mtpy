import matplotlib.pyplot as plt
import numpy as np


def get_next_fig_num():
    current_fig_nums = set(plt.get_fignums())
    number = 1
    while number in current_fig_nums:
        number += 1
    return number


def gen_hist_bins(uniq_period_list):
    bins = np.array(uniq_period_list)  # get center of bins
    diff = np.diff(np.r_[0, bins]).dot(.5)
    bins -= diff  # shift left
    bins = np.r_[bins, uniq_period_list[-1] + diff[-1]]  # add last bar
    return bins
