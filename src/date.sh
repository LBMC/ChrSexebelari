#!/bin/sh

# This script renames files in a given directory using git repository file_handle. 
# bash date.sh "dir"

PATH_FILES=$1
FILES="$(ls $PATH_FILES)"

for f in ${FILES} 
do
	python3 ./src/file_handle/src/file_handle.py -f $PATH_FILES/$f
done
