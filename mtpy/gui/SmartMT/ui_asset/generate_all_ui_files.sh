#!/usr/bin/env bash

#ls -1 *.ui | while read file
for afile in `ls *.ui`;
do
    afilebn=`basename $afile .ui`
    echo "generating $afilebn.py ......"
    pyuic4 $afile > $afilebn.py
done
