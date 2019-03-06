from . import qel_prepare_birrp_data3 as ppb
import os
import sys
import os.path as op
import time

station_names = ['01']

# station_names = ['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19',
# '20','21','22','23','24','25']


lo_subfolders = ['/media/Elements2/KEX01_32gig/']


dates = ['2014-01-30-at-13-00-00']
# ,'2014-01-17-at-13-00-00','2014-01-18-at-13-00-00',
# '2014-01-19-at-13-00-00','2014-01-20-at-13-00-00','2014-01-21-at-13-00-00','2014-01-22-at-13-00-00',
# '2014-01-23-at-13-00-00','2014-01-24-at-13-00-00','2014-01-25-at-13-00-00','2014-01-26-at-13-00-00','2014-01-27-at-13-00-00',
# '2014-01-28-at-13-00-00','2014-01-29-at-13-00-00','2014-01-30-at-13-00-00','2014-01-31-at-13-00-00','2014-02-01-at-13-00-00',
# '2014-02-02-at-13-00-00','2014-02-03-at-13-00-00','2014-02-04-at-13-00-00','2014-02-05-at-13-00-00','2014-02-06-at-13-00-00',
# '2014-02-07-at-13-00-00','2014-02-08-at-13-00-00','2014-02-09-at-13-00-00','2014-02-10-at-13-00-00',
# '2014-02-11-at-13-00-00','2014-02-12-at-13-00-00','2014-02-13-at-13-00-00','2014-02-14-at-13-00-00', '2014-02-15-at-13-00-00',
# '2014-02-16-at-13-00-00','2014-02-17-at-13-00-00','2014-02-18-at-13-00-00',
# '2014-02-19-at-13-00-00','2014-02-20-at-13-00-00','2014-02-21-at-13-00-00','2014-02-22-at-13-00-00','2014-02-23-at-13-00-00',
# '2014-02-24-at-13-00-00','2014-02-25-at-13-00-00','2014-02-26-at-13-00-00','2014-02-27-at-13-00-00','2014-02-28-at-13-00-00',
# '2014-03-01-at-13-00-00','2014-03-02-at-13-00-00','2014-03-03-at-13-00-00','2014-03-04-at-13-00-00','2014-03-05-at-13-00-00',
# '2014-03-06-at-13-00-00','2014-03-07-at-13-00-00','2014-03-08-at-13-00-00',
# '2014-03-09-at-13-00-00','2014-03-10-at-13-00-00','2014-03-11-at-13-00-00','2014-03-12-at-13-00-00','2014-03-13-at-13-00-00',
# '2014-03-14-at-13-00-00','2014-03-15-at-13-00-00']


channels = [(0, 1)]
# channels = [(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),
# (0,1),(0,1),(0,1),(0,1),(0,1),(0,1)]


outputbase = '/stash/elogger/L101_30_Jan_Test/400000_detrend_removed'

channelsB = (0, 1)

dataDirB = '/media/Elements/KEXRR_32gig/RRB/'

channelsC = (0, 1)

dataDirC = '/media/Elements2/KEX01_32gig/L125_Blogger/'

prefix = 'L1'

subfolder_addon = '_RR_B125'

for idx_s, station in enumerate(station_names):

    dataDirA = op.join(lo_subfolders[0], prefix + station)

    subfoldername = '{0}{1}{2}'.format(prefix, station, subfolder_addon)

    subfolder = op.join(outputbase, subfoldername)
    if not op.isdir(subfolder):
        os.makedirs(subfolder)

    channelsA = channels[idx_s]

    print(dataDirA, dataDirB, dataDirC)
    print(channelsA, channelsB, channelsC)

    for idx_date, date in enumerate(dates):
        timestamp = time.strptime(date, "%Y-%m-%d-at-%H-%M-%S")
        try:
            # print 	date,subfolder
            ppb.run(
                dataDirA,
                dataDirB,
                dataDirC,
                channelsA,
                channelsB,
                channelsC,
                timestamp,
                subfolder)
            # sys.exit()
        except:
            continue
