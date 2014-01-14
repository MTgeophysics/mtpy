#!/usr/bin/env python



import numpy as np
from pylab import *



def merge_periods(periods,merge_threshold):
    """
    assume periods in increasing order
    merge_threshold given in percent

    """
    old_periods = sorted(list(periods),reverse=False)

    in_values = [(x) for x in old_periods]


    out_values = []
    cluster_counter = 0
    periods_in_cluster = []
    finished_cluster = False
    cluster_flags = np.zeros_like(in_values)+1 
    
    idx = 0 
    while idx < len(in_values):


        base_period = in_values[idx]

        while finished_cluster is False:
            

            periods_in_cluster.append(base_period)
            if idx != 0 :
                cluster_flags[idx] = cluster_counter
    
            if idx == len(in_values) - 1:
                finished_cluster = True
                continue


            #check, if the next period is too far away from the current mean:           
            mean_period = np.mean(periods_in_cluster)
            distance = in_values[idx+1] - mean_period


            if distance/mean_period >  merge_threshold/100.:
                #print 'next point too far...finishing cluster'
                finished_cluster = True
                continue
        
            #check, if the first period would drop out of the threshold range,
            # if the cluster gets extended:
            extended_periods = periods_in_cluster[:]
            extended_periods.append(in_values[idx+1])
            new_mean = np.mean(extended_periods) 
            distance_period1 = new_mean - extended_periods[0]


            if distance_period1/new_mean > merge_threshold/100.:
                #print 'first point would drop out...finishing cluster'
                finished_cluster = True
                continue
            #otherwise just continue filling the same cluster
            idx += 1

        if finished_cluster is True:
            mean_period = [np.mean(periods_in_cluster)]
            if len(periods_in_cluster) > 1:
                #print 'merging periods', periods_in_cluster,'into' , mean_period 
                pass

            out_values.extend(len(periods_in_cluster)*mean_period)
            finished_cluster = False
            periods_in_cluster = []
            idx += 1
            cluster_counter += 1


    new_period_list = out_values
    if len(in_values) != cluster_counter:
        print '\n\tDone -- merged {0} periods into {1} period-clusters\n'.format(len(in_values),cluster_counter)
    return new_period_list


def plot_merging(periods,merge_threshold):


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

    mergedperiods = merge_periods(periods,merge_threshold)
    



    ax = subplot2grid((1, 1), (0, 0), colspan=1)
    orig = ax.scatter(periods,zeros(len(periods)),label='original periods')
    #ax.set_xlim([10**(-rng),10**(rng)])
    ax.set_xscale('log')

    hold(True)

    mergedperiods = sorted(list(set(mergedperiods)),reverse=False)

    lo_limits = []
    for i,p in enumerate(mergedperiods):
        if i == 0:
            continue
        lo_limits.append(np.sqrt(mergedperiods[i] *mergedperiods[i-1]))


    ax.set_xticks(lo_limits, minor=True)
    ax.xaxis.grid(False, which='major')
    ax.xaxis.grid(True, which='minor',c='g')

    merge = ax.scatter(mergedperiods,ones(len(mergedperiods)),c='r',label='merged periods')
    ax.set_ylim([-1,2])
    handles, labels = ax.get_legend_handles_labels()
    ax.legend([orig, merge], ["original ({0})".format(len(periods)), "merged ({0})".format(len(mergedperiods))],scatterpoints=1,loc='upper center', 
                ncol=2)
    
    #tight_layout()
    show()#block=True)
    raw_input()



if __name__ == '__main__':


    N = 100
    rng= 4
    rnd = np.array(sorted((np.random.random_sample(N)-0.5)*2*rng,reverse=False))
    periods = 10**(rnd)
    threshold = 5
    print """

        This is a module - not to be run as as a script! 

This call yields an example result plot for merging {0} random periods.

""".format(N)

    plot_merging(periods,threshold)

