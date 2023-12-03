import os,sys
import pysam

def long_short_split(sample_path,sample_name_list):
    output_path = os.path.join(sample_path,'fragment_study/mis_10_splitedBy125','125_splited_sam_file')
#    problem_file_path = os.path.join(sample_path,'mis_10_splitedBy125','problem_file')
#     os.chdir(sample_path)
    os.makedirs(output_path,exist_ok=True)
#    os.makdirs(problem_file_path)
    sample_name_list_path = os.path.join(sample_path,sample_name_list)
    f = open(sample_name_list_path)
    sample_names = f.readlines()
    total_sample = 0
    for i in sample_names:
        i = i.replace('\n','').strip()
        if os.path.getsize(os.path.join(sample_path,i,i+'.mis.10.bam')) > 1:
            item_sam_file_path = os.path.join(sample_path,i,i+'.mis.10.bam')
            print('Processing file '+ i+'.mis.10.bam')
            sam_fitered_by_insertSize_path4 = os.path.join(output_path,i+'.mis.10.S125.sam')
            sam_fitered_by_insertSize_path5 = os.path.join(output_path,i+'.mis.10.L125.sam')
            sam_ali = pysam.AlignmentFile(item_sam_file_path,'rb')

            sam_fitered_by_insertSize4 = pysam.AlignmentFile(sam_fitered_by_insertSize_path4,'w',template = sam_ali)
            sam_fitered_by_insertSize5 = pysam.AlignmentFile(sam_fitered_by_insertSize_path5,'w',template = sam_ali)
            sam_ali_fetch = sam_ali.fetch()
            for read in sam_ali_fetch:
                insert_size = abs(int(read.template_length))
                if insert_size<125:#urine135,plasma125
                    sam_fitered_by_insertSize4.write(read)
                elif insert_size>=125:
                    sam_fitered_by_insertSize5.write(read)
            sam_ali.close()
            sam_fitered_by_insertSize4.close()
            sam_fitered_by_insertSize5.close()
            total_sample += 1
            print(str(total_sample)+'samples have been processed')
            if total_sample == len(sample_names):
                print('All samples have been processed')
        else:
#            pf = open(os.path.join(problem_file_path,'sample_no_mis10.txt','a+')
#            pf.write(i+'\n')
            print('File '+i+' mis10 too small')
        #print("cal depth!")
        #os.system(f"samtools depth -d 500000000 -a -m 0 ${sam_fitered_by_insertSize_path4} > ../depth_file/${s}.txt"
            
    f.close()
           
def run():
    sample_path = sys.argv[1]
    #sample_path = "/mnt/data1/Seq_data/20230815/N_Plasma_DOU_mtDNA-CAP/Analysis_private-mtDNA_V1"
    sample_name_list = 'sample_name.txt'
    long_short_split(sample_path,sample_name_list)
if __name__ == '__main__':
    run()
