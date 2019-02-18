#!/usr/bin/env python


from numpy import *
import sys
import os
import os.path as op

channel_dict = {'n': 0, 'e': 1, 's': 2, 'w': 3}
longchannel_dict = {'north': 0, 'east': 1, 'south': 2, 'west': 3,
                    '0': 'north', '1': 'east', '2': 'south', '3': 'west'}


def main():

    if len(sys.argv) < 7:
        sys.exit('\n\tERROR - need 6 arguments as input: \n <input dir> <output dir>'
                 ' <stationA> <stationB> <stationC> <prefix>\n')
    print()

    indir = sys.argv[1]
    indir = op.abspath(op.join(os.curdir, indir))

    if not op.isdir(indir):
        sys.exit('\n\t ERROR - input directory does not exist! \n')

    outdir = sys.argv[2]
    outdir = op.abspath(op.join(os.curdir, outdir))

    stationnames = [i.upper() for i in sys.argv[3:6]]

    prefix = sys.argv[6]

    t0, sampling = read_timestampfile(indir, prefix)

    if not op.isdir(outdir):
        os.makedirs(outdir)

    lo_infiles = os.listdir(indir)
    lo_infiles = [
        i for i in lo_infiles if i.lower().startswith(
            prefix.lower())]

    for idx_sta, sta in enumerate(['A', 'B', 'C']):
        fullprefix = prefix + '.sta' + sta
        station = stationnames[idx_sta]
        outfn = op.join(outdir, '%s_1Hzoverview_fd' % station)
        baseline = []
        alldata = []
        for c in range(4):
            infile = op.join(indir, fullprefix + '.%s' %
                             (longchannel_dict[str(c)]))
            tmp_data = []
            Fin = open(infile)
            firstvalue = int(float(Fin.readline().strip().split()[0]))

            baseline.append(firstvalue)
            tmp_data.append(0)
            oldval = firstvalue
            for line in Fin.readlines():
                val = int(float(line.strip().split()[0]))
                tmp_data.append(val - oldval)
                oldval = val
            alldata.append(tmp_data)
            Fin.close()
        Fout = open(outfn, 'w')
        header = '# {0:d} \t{1:d} \t{2:d} \t{3:d} \t\t{4:.4f} \n'.format(baseline[0],
                                                                         baseline[1], baseline[2], baseline[3], t0)
        Fout.write(header)
        for i in range(len(alldata[0])):
            s = '{0:d} \t{1:d} \t{2:d} \t{3:d}\n'.format(alldata[0][i], alldata[1][i],
                                                         alldata[2][i], alldata[3][i])
            Fout.write(s)
        Fout.close()
        print('\tdata file written: {0}\n'.format(outfn))

    print()
    print('\tDone!\n')


def read_timestampfile(indir, prefix):

    fn = prefix + '.timestamps'
    fn = op.join(indir, fn)
    t0 = None
    sampling = None
    in_dict = {}
    Fin = open(fn)
    for line in Fin:
        line = line.strip()
        try:
            in_dict[line.split(':')[0].strip()] = line.split(':')[1].strip()
        except:
            continue

    for k, v in list(in_dict.items()):
        if k.lower().startswith('sampl'):
            sampling = float(v)
            if sampling % 1 == 0:
                sampling = int(sampling)

        if k.lower().startswith('first'):
            t0 = float64(v)

    return t0, sampling


if __name__ == '__main__':
    main()
