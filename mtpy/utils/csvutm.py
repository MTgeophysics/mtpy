#!/usr/bin/env python
'''Convert and add columns for different coordinate systems to a CSV file.

This script requires pyproj installed. If you have a CSV file with two columns
containing x and y coordinates in some coordinate system (the "from" system),
this script will add another two columns with the coordinates transformed into
another system (the "to" system). The CSV file must have a header row, and
the column names should be specified on the command line as well.

The coordinate systems are defined by their EPSG codes. Geographic coordinates
(decimal degrees) on the WGS-84 datum are code 4326. There is a code for
every system imaginable: see a searchable list at

http://spatialreference.org/ref/epsg/

'''
# Standard library packages
import argparse
import sys
try:
    import cStringIO as StringIO
except ImportError:
    import StringIO

# Third-party packages
import csv
import pyproj


def main():
    parser = get_parser()
    args = parser.parse_args(sys.argv[1:])
    with open(args.in_csv_filename[0], mode='r') as f:
        csv_txt = f.read()
    out_file = open(args.out_csv_filename[0], mode='wb')
    csvutm(csv_txt, out_file, delimiter=args.delimiter,
           f=args.from_coords, fx=args.fx, fy=args.fy,
           t=args.to, tx=args.tx, ty=args.ty)
    

def csvutm(csvtxt, out_file, delimiter=',',
           f='28353', fx='easting', fy='northing',
           t='4326', tx='lon', ty='lon'):

    """
    ...

    """
    f_in = StringIO.StringIO(csvtxt)
    r = csv.DictReader(f_in, delimiter=delimiter)
    for key in (fx, fy):
        try:
            assert key in r.fieldnames
        except AssertionError:
            print('Did not find %s in CSV file.' % key)
            raise
    fxs = []
    fys = []
    for i, row in enumerate(r):
        fxs.append(float(row[fx]))
        fys.append(float(row[fy]))
    p1 = pyproj.Proj(init='epsg:%s' % f)
    p2 = pyproj.Proj(init='epsg:%s' % t)
    txs, tys = pyproj.transform(p1, p2, fxs, fys)
    f_in.seek(0)
    r = csv.DictReader(f_in, delimiter=delimiter)
    fnames = r.fieldnames
    if tx not in fnames:
        fnames += [tx]
    if ty not in fnames:
        fnames += [ty]
    w = csv.DictWriter(out_file, fieldnames=fnames, delimiter=delimiter)
    w.writeheader()
    for i, row in enumerate(r):
        row[tx] = txs[i]
        row[ty] = tys[i]
        w.writerow(row)

            
def get_parser():
    """
    ...

    """
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0],
            epilog=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--fx', default='lon', help='column header for x coord of 1st (from) coord system')
    parser.add_argument('--fy', default='lat', help='column header for y coord of 1st (from) coord system')
    parser.add_argument('--tx', default='easting', help='column header for x coord of 2nd (to) coord system')
    parser.add_argument('--ty', default='northing', help='column header for y coord of 2nd (to) coord system')
    parser.add_argument('-f', '--from', help='EPSG code for coordinate system to convert from.\n'
                                             'See http://spatialreference.org/ref/epsg/', default='4326', dest='from_coords')
    parser.add_argument('-t', '--to', help='EPSG code for coordinate system to convert into.', default='28353')
    parser.add_argument('-d', '--delimiter', default=',')
    parser.add_argument('in_csv_filename', nargs=1)
    parser.add_argument('out_csv_filename', nargs=1)
    return parser
    
    
if __name__ == '__main__':
    main()
