sample_dir=$1
cd ${sample_dir}/fragment_study/mis_10_splitedBy125
mkdir depth_file
cd ${sample_dir}/fragment_study/mis_10_splitedBy125/125_splited_sam_file
fils=$(ls -l *.sam | awk '{if($5>1000) print$9}')
for s in ${fils}
do
    samtools depth -d 500000000 -a -m 0 ${s} > ../depth_file/${s}.txt
done