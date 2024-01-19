#!/bin/bash

#rootfiles="/data/11b/itabashi/root/*"
rootfiles="../current_info/*"
rm -rf current_files.list

for pathfile in $rootfiles; do
	echo $pathfile >> "current_files.list"
done
	
