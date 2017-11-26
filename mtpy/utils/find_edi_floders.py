# -*- coding: utf-8 -*-
# ! /usr/bin/env python
"""
Description:
    Find path to directories which containing a given type of files: .edi, .py .jpg, .pdf


CreationDate:   26/11/2017
Developer:      fei.zhang@ga.gov.au

Revision History:
    LastUpdate:  26/11/2017   FZ started the first version
"""

import os
import sys

class Search4EdiFolders(object):

    def __init__(self, startDir, edifiles_threshold=0, filetype='edi'):
        self.startDir = startDir
        self.edifiles_threshold =edifiles_threshold
        self.filetype =filetype
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
                if fname.lower().endswith(self.filetype):
                #if fname.endswith('.py'):
                    edi_files_count = edi_files_count+1
                else:
                    pass

            if edi_files_count > self.edifiles_threshold:
                print('Found directory: %s ==> %s *.%s files' %(dirName, edi_files_count,self.filetype))
                self.edi_folders.append(dirName)

            return self.edi_folders


########################################################, edi_files_count
if __name__ =="__main__":

    if len(sys.argv)<=1:
        print("USAGE: %s  %s [%s] [%s]" % (sys.argv[0], "start_path_dir", "ftype", "edi_threshold"))
        sys.exit(1)
    else:
        root_dir=sys.argv[1]

        if len(sys.argv)>2:
            ftype = sys.argv[2].lower()
        else:
            ftype='edi'

        if len(sys.argv)>3:
            edifiles_th = int(sys.argv[3])
        else:
            edifiles_th=0


        sObj= Search4EdiFolders(root_dir, edifiles_threshold=edifiles_th, filetype=ftype)
        #sObj= Search4EdiFolders(root_dir, edifiles_threshold=edifiles_th, filetype='py')

        edi_dirs_list =sObj.find_edi_folders(root_dir)

        print (edi_dirs_list)
        print ("Number of interesting folders found =", len(edi_dirs_list))
