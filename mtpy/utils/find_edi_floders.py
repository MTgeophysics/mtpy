# -*- coding: utf-8 -*-
# ! /usr/bin/env python
"""
Description:
    Find EDI files-containing directories - their full path

CreationDate:   26/11/2017
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:  26/11/2017   FZ started the first version
"""

import os
import sys

class Search4EdiFolders(object):

    def __init__(self, startDir, edifiles_threshold=0):
        self.startDir = startDir
        self.edifiles_threshold =edifiles_threshold
        self.edi_folders=[] # initial empty list

    def find_edi_folders(self, aStartDir):
        """ find edi-files containing directories
        rootDir: the directory you want to start from
        """
        for dirName, subdirList, fileList in os.walk(aStartDir):

            for asubdir in subdirList:
                self.find_edi_folders(os.path.join(dirName,asubdir))

            edi_files_count=0
            for fname in fileList:
                if fname.endswith('edi'):
                    edi_files_count = edi_files_count+1
                else:
                    pass

            if edi_files_count > self.edifiles_threshold:
                print('Found EDI directory: %s ==> %s EDI files' %(dirName, edi_files_count))
                self.edi_folders.append(dirName)

            return self.edi_folders


########################################################
if __name__ =="__main__":

    if len(sys.argv)<=1:
        print("USAGE: %s  %s [%s]" % (sys.argv[0], "start_path_dir", "edi_threshold"))
        sys.exit(1)
    else:
        root_dir=sys.argv[1]

        if len(sys.argv)>2:
            edifiles_th = int(sys.argv[2])
        else:
            edifiles_th=0


        sObj= Search4EdiFolders(root_dir, edifiles_threshold=edifiles_th)

        edi_dirs_list =sObj.find_edi_folders(root_dir)

        print (edi_dirs_list)
        print ("Number of EDI_folders found =", len(edi_dirs_list))
