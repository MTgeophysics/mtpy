# -*- coding: utf-8 -*-
# ! /usr/bin/env python
"""
Description:
    Find path to all the directories which contain a given type of files: .edi, .py .jpg, .pdf

How to Run:
    python mtpy/utils/edi_folders.py  . EDI
    python mtpy/utils/edi_folders.py  /e/Data/ EDI 2
    python mtpy/utils/edi_folders.py  /e/Data/ PY

CreationDate:   26/11/2017
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:  26/11/2017   FZ started the first version

"""

import os
import sys


def recursive_glob(dirname, ext='*.edi'):
    """
    Under the dirname recursively find all files with extension ext.
    Return a list of the full-path to the types of files of interest.

    This function is useful to handle a nested directories of EDI files.

    :param dirname: a single dir OR a list of dirs.
    :param ext: eg, ".edi", ".xml"
    :return: a list of path2files
    """
    import fnmatch

    if isinstance(dirname, list):  # the input argument is a list of directory
        filelist = []
        for adir in dirname:
            filelist.extend(recursive_glob(adir))
        return filelist
    else:  # input variable is a single dir
        matches = []
        for root, dirnames, filenames in os.walk(dirname):
            for filename in fnmatch.filter(filenames, ext):
                matches.append(os.path.join(root, filename))
        return matches


class EdiFolders(object):

    def __init__(self, startDir, edifiles_threshold=1, filetype='.edi'):
        self.startDir = startDir  # the top level dir to be searched
        self.edifiles_threshold = edifiles_threshold  # at least 1 file is of the specified type.
        self.filetype = filetype
        self.folders_of_interest = []  # initial empty list

    def find_edi_folders(self, aStartDir):
        """
        find  directories containing the file of type self.filetype
        :param aStartDir: the directory to start from
        :return: a list of full path to folders of interest.
        """

        for dirName, subdirList, fileList in os.walk(aStartDir):

            edi_files_count = 0
            for fname in fileList:
                if fname.lower().endswith(self.filetype):
                    # if fname.endswith('.py'):
                    edi_files_count = edi_files_count + 1
                else:
                    pass

            if edi_files_count >= self.edifiles_threshold:
                print(('Found directory: %s ==> %s *.%s files' % (dirName, edi_files_count, self.filetype)))
                self.folders_of_interest.append(dirName)

            #If it's a folder then recursive call
            for asubdir in subdirList:
                self.find_edi_folders(os.path.join(dirName, asubdir))

            return self.folders_of_interest

    def get_all_edi_files(self):

        return recursive_glob(self.startDir)


########################################################
if __name__ == "__main__":

    if len(sys.argv) <= 1:
        print(("USAGE: %s  %s [%s] [%s]" % (sys.argv[0], "start_path_dir", "ftype", "edi_threshold")))
        sys.exit(1)
    else:
        root_dir = sys.argv[1]

        if len(sys.argv) > 2:
            ftype = sys.argv[2].lower()
        else:
            ftype = 'edi'

        if len(sys.argv) > 3:
            edifiles_th = int(sys.argv[3])
        else:
            edifiles_th = 1

        sObj = EdiFolders(root_dir, edifiles_threshold=edifiles_th, filetype=ftype)

        edi_dirs_list = sObj.find_edi_folders(root_dir)

        print(("Number of interesting folders found =", len(edi_dirs_list)))

        print (edi_dirs_list)
