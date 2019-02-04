#!/usr/bin/env python

"""
    remove_bz_files.py

    Script for removing all files ending in *.BZ in a directory (including subfolders)
    Useful, if BZ channel data is just junk output of the data logger, but it's
    included in the general filehandling by default.

"""

import os
import sys
import os.path as op


def main():

    if len(sys.argv) < 2:
        print('\nNeed at least 1 argument: <directory to clean> \n\n'\
            'Optional flag: \n [-R]\n'\
            ' (recursive)\n\n')
        return

    path = sys.argv[1]
    try:
        path = op.abspath(op.join(os.curdir, path))
        if not op.isdir(path):
            raise
    except:
        sys.exit('Data file(s) path not existing: {0}\n'.format(path))

    rec = False

    if len(sys.argv) > 2:
        optionals = sys.argv[2:]
        for o in optionals:
            o = o.replace('-', '')
            if o.lower().startswith('r'):
                rec = True

    remove_files(path, rec)


def remove_files(path, recursive_flag=False):

    if recursive_flag is not True:
        lo_files = os.listdir(path)
        lo_files = [op.join(path, i) for i in lo_files]
        lo_files = [i for i in lo_files if op.isfile(i)]
        lo_files = [i for i in lo_files if i.lower().endswith('.bz')]

    else:
        lo_files = [op.join(dp, f) for dp, dn, filenames in os.walk(
            path) for f in filenames if op.splitext(f.lower())[1] == '.bz']

    if len(lo_files) == 0:
        print('\nFound no files to delete\n')
        return

    print('\nFound files to delete:\n {0}'.format(lo_files))

    confirm = False

    while confirm is False:
        answer = input(
            '\n\t Do you want to remove the files permanently? (y/n) ')
        try:
            answer = answer.lower()[0]
            if answer in ['y', 'n']:
                confirm = True
        except:
            pass

    print()
    if answer == 'y':
        print('....deleting files...', end=' ')
        for f in lo_files:
            os.remove(f)
        print('Done!\n')

    return


if __name__ == '__main__':
    main()
