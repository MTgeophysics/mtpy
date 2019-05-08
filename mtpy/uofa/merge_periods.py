#!/usr/bin/env python


import numpy as np
from pylab import *


def regular_periods(periodlist, merge_threshold=15, no_periods=None, t_min=None,
                    t_max=None, max_merge_error=None):
    """merging bins are rgegularly shaped on the log10 axis around the respective center frequencies

    """

    if no_periods is None:
        no_periods = 20
    if t_min is None or min(periodlist) > t_min:
        t_min = min(periodlist)
    if t_max is None or max(periodlist) < t_max:
        t_max = max(periodlist)
    if max_merge_error is None:
        max_merge_error = 10.

    new_periods = np.logspace(np.log10(t_min), np.log10(t_max), no_periods)

    new_periods_log = log10(new_periods)
    bin_width = merge_threshold / 100. * \
        (max(new_periods_log) - min(new_periods_log))
    new_periods = np.logspace(
        np.log10(t_min) +
        0.25 *
        bin_width,
        np.log10(t_max) -
        0.25 *
        bin_width,
        no_periods)

    new_period_bins = []

    for log_p in new_periods_log:
        min_p = 10**(log_p - bin_width / 2.)
        max_p = 10**(log_p + bin_width / 2.)
        new_period_bins.append([min_p, max_p])

    new_period_list = []
    merge_errors = []

    for p in periodlist:
        new_p = None
        p_error = None
        for idx_p, p_bin in enumerate(new_period_bins):
            if p_bin[0] <= p <= p_bin[1]:
                new_p = new_periods[idx_p]
                deviation = np.abs(p - new_p)
                try:
                    max_deviation = max(
                        np.abs(
                            new_p -
                            p_bin[0]),
                        np.abs(
                            p_bin[1]) -
                        new_p)
                    p_error = deviation / max_deviation
                except:
                    # only possible, if bin-width is zero, so it's exactly
                    # correct:
                    p_error = 0.

        new_period_list.append(new_p)
        merge_errors.append(p_error)

    merged_periods = set([i for i in new_period_list if i is not None])

    ignored_points = [True for i in new_period_list if i is None]

    if len(merged_periods) != len(periodlist):
        print('\n\tMerged {0} periods into {1} period-clusters -'\
            ' {2} points outside the bins\n'.format(
                len(periodlist), len(merged_periods), len(ignored_points)))

    new_period_list = [
        round(
            i,
            5) if i is not None else i for i in new_period_list]

    return new_period_list, merge_errors


def merge_periods(periods, merge_threshold):
    """
    - assume that the periods are in increasing order
    - merge_threshold given in percent

    """
    old_periods = sorted(list(periods), reverse=False)

    in_values = [(x) for x in old_periods]

    out_values = []
    cluster_counter = 0
    periods_in_cluster = []
    finished_cluster = False
    cluster_flags = np.zeros_like(in_values) + 1

    idx = 0
    while idx < len(in_values):

        base_period = in_values[idx]

        while finished_cluster is False:

            periods_in_cluster.append(base_period)
            if idx != 0:
                cluster_flags[idx] = cluster_counter

            if idx == len(in_values) - 1:
                finished_cluster = True
                continue

            # check, if the next period is too far away from the current mean:
            mean_period = np.mean(periods_in_cluster)
            distance = in_values[idx + 1] - mean_period

            if distance / mean_period > merge_threshold / 100.:
                # print 'next point too far...finishing cluster'
                finished_cluster = True
                continue

            # check, if the first period would drop out of the threshold range,
            # if the cluster gets extended:
            extended_periods = periods_in_cluster[:]
            extended_periods.append(in_values[idx + 1])
            new_mean = np.mean(extended_periods)
            distance_period1 = new_mean - extended_periods[0]

            if distance_period1 / new_mean > merge_threshold / 100.:
                # print 'first point would drop out...finishing cluster'
                finished_cluster = True
                continue
            # otherwise just continue filling the same cluster
            idx += 1

        if finished_cluster is True:
            mean_period = [np.mean(periods_in_cluster)]
            if len(periods_in_cluster) > 1:
                # print 'merging periods', periods_in_cluster,'into' ,
                # mean_period
                pass

            out_values.extend(len(periods_in_cluster) * mean_period)
            finished_cluster = False
            periods_in_cluster = []
            idx += 1
            cluster_counter += 1

    new_period_list = out_values
    if len(in_values) != cluster_counter:
        print('\n\tDone -- merged {0} periods into {1} period-clusters\n'.format(len(in_values), cluster_counter))
    return new_period_list


def plot_merging(periods, merge_threshold, no_periods=None, t_min=None,
                 t_max=None, max_merge_error=None):

    # import platform,os,sys
    # if not platform.system().lower().startswith('win') :

    #     #generate an interactive plot window, which remains open after this script has finshed:
    #     proc_num = os.fork()

    #     if proc_num != 0:
    #         #This is the parent process, that should quit immediately to return to the
    #         #shell.
    #         print "You can kill the plot window with the command \"kill %d\"." % proc_num
    #         sys.exit()

    close('all')
    ion()

    if no_periods is not None:
        mergedperiods, dummy = regular_periods(periods, merge_threshold, no_periods,
                                               t_min, t_max, max_merge_error)
    else:
        mergedperiods = merge_periods(periods, merge_threshold)

    ax = subplot2grid((1, 1), (0, 0), colspan=1)
    orig = ax.scatter(periods, zeros(len(periods)), label='original periods')
    # ax.set_xlim([10**(-rng),10**(rng)])
    ax.set_xscale('log', nonposx='clip')

    hold(True)

    n_points_used = len([i for i in mergedperiods if i is not None])
    n_points_outside = len([i for i in mergedperiods if i is None])

    mergedperiods = [i for i in mergedperiods if i is not None]

    mergedperiods = sorted(list(set(mergedperiods)), reverse=False)

    lo_limits = []

    new_periods_log = log10(np.array(mergedperiods))
    bin_width = merge_threshold / 100. * \
        (max(new_periods_log) - min(new_periods_log))

    new_period_bins = []

    for log_p in new_periods_log:
        min_p = 10**(log_p - bin_width / 2.)
        max_p = 10**(log_p + bin_width / 2.)
        lo_limits.append(min_p)
        lo_limits.append(max_p)

    ax.set_xticks(lo_limits, minor=True)
    ax.xaxis.grid(False, which='major')
    ax.xaxis.grid(True, which='minor', c='g')

    merge = ax.scatter(
        mergedperiods,
        ones(
            len(mergedperiods)),
        c='r',
        label='merged periods')
    ax.set_ylim([-1, 2])
    ax.set_xlim([10**(min(new_periods_log) - 0.5),
                 10**(max(new_periods_log) + 0.5)])
    handles, labels = ax.get_legend_handles_labels()
    ax.legend([orig, merge], ["original ({0})".format(len(periods)), "merged ({0})".format(len(mergedperiods))], scatterpoints=1, loc='upper center',
              ncol=2)
    ax.set_title(
        '{0} periods in bins - {1} periods left out'.format(n_points_used, n_points_outside))

    tight_layout()
    show()  # block=True)
    input()


if __name__ == '__main__':

    N = 100
    rng = 4
    rnd = np.array(
        sorted(
            (np.random.random_sample(N) - 0.5) * 2 * rng,
            reverse=False))
    periods = 10**(rnd)
    threshold = 10
    print("""

        This is a module - not to be run as as a script!

This call yields an example result plot for merging {0} random periods.

""".format(N))

    no_periods = 5
    print(min(periods), max(periods))

    plot_merging(periods, threshold, no_periods)
