#!bin/bash
start=$(date +%s) #start time
echo "start time: "`date` 
sed -i 's/\r//g' sample_path_list.txt
cat sample_path_list.txt | while read line
do

    inpath=${line}
    echo ${inpath}
    #cal endmotifs
    #python /mnt/hdd/usr/jiaohuanmin/home/jiaohm/Tools/mtDNA_pipeline/fragment_stduy_pipeline/endmotifs_stat_v1.py ${inpath}
    # cal LS125 sam
    #python /mnt/hdd/usr/jiaohuanmin/home/jiaohm/Tools/mtDNA_pipeline/fragment_stduy_pipeline/split_mis10_by125.py ${inpath}
    #cal depth
    #bash /mnt/hdd/usr/jiaohuanmin/home/jiaohm/Tools/mtDNA_pipeline/fragment_stduy_pipeline/chrM_depth.sh ${inpath}
    #cal zfsd
    #python /mnt/hdd/usr/jiaohuanmin/home/jiaohm/Tools/mtDNA_pipeline/fragment_stduy_pipeline/concate_depth_v4.py ${inpath}
    #cal zfsd feature
    #python /mnt/hdd/usr/jiaohuanmin/home/jiaohm/Tools/mtDNA_pipeline/fragment_stduy_pipeline/fsd_feature_cal_v0.1.py ${inpath}
    # cal insertsize
    #bash /mnt/hdd/usr/jiaohuanmin/home/jiaohm/Tools/mtDNA_pipeline/fragment_stduy_pipeline/insertisize_find_mt.sh ${inpath}
    #concate insertsize
    python /mnt/hdd/usr/jiaohuanmin/home/jiaohm/Tools/mtDNA_pipeline/fragment_stduy_pipeline/insertsize_concate_v3.1.py ${inpath}
done
