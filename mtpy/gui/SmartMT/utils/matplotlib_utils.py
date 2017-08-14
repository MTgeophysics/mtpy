import matplotlib.pyplot as plt


def get_next_fig_num():
    current_fig_nums = set(plt.get_fignums())
    number = 1
    while number in current_fig_nums:
        number += 1
    return number
