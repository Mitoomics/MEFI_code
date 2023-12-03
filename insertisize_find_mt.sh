#!bin/bash
result_dir=$1
#result_dir=/mnt/data1/jiaohm/mapping_methods/Analysis_rCRS-hg38-DNA_V1
cd $result_dir
mkdir tmp
filenames=$(ls ./*/*mis.*.bam)
for i in ${filenames}
do
    name=${i##*/}
    sample=${name%%.*}
    echo ${name}
    java -Xmx10g -Djava.io.tmpdir=`pwd`/tmp -jar /mnt/hdd/softs/picard/picard-tools-1.81/CollectInsertSizeMetrics.jar I=${i} O=./${sample}/${name}.insertsize.txt H=./${sample}/${name}.insertsize.pdf VALIDATION_STRINGENCY=LENIENT 2>>./insertsize.analysis.log
done
rm -rf tmp
