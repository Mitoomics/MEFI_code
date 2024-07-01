#!/bin/bash
start=$(date +%s) #start time
echo "start time: "`date` 

if [ ! -f sample_path_list.txt ]; then
   echo "sample_path_list.txt not exist, check again!"
   exit 1
fi

sed -i 's/\r//g' sample_path_list.txt

cat sample_path_list.txt | while read line
do
    echo ${line}
    #cal endmotifs
    #python endmotifs_stat.py ${line}
    # cal LS125 sam
    #python split_mis10_by125.py ${line}
    #cal depth
    #bash chrM_depth.sh ${line}
    #cal fsd and fsd feature
    #python fsd_compute.py ${line}
    # compute insertsize
    #bash insertisize_find_mt.sh ${line}
    #insertsize stat
    python insertsize_stat.py ${line}
done
