#!/usr/bin/env bash

ls -1 *.ui | while read file
do
    set x
    pyuic '$file' > '$(basename '$file').py'
done