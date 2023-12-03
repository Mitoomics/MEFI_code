import os, sys
import glob
import pandas as pd
from median_cal import median_isnertsize


class FileConcate:
    def __init__(self,
                 sample_path,
                 file_name,
                 target_path=None,
                 specific_path=None):
        self.sample_path = sample_path
        self.target_path = target_path
        self.specific_path = specific_path
        self.file_name = file_name
        #self.result_path = result_path

    def create_insertsize_mtr(self):
        df = pd.DataFrame({'insert_size': range(1001)})
        df.to_csv(self.result_file, index=False)
        print(f'created {file_name} matrix done!')

    def file_concate(self):
        outpath = os.path.join(self.sample_path, "fragment_study/insertszie")
        os.makedirs(outpath, exist_ok=True)
        self.result_file = os.path.join(outpath,
                                        self.file_name)  #save to sample root
        #self.result_file = os.path.join(i,result_path,self.file_name) # save to destinational root
        self.create_insertsize_mtr()
        sam_path = os.path.join(self.sample_path,
                                '*/*.mis.10.bam.insertsize.txt')
        #                 specific_file_path =  os.path.join(i,self.target_path,self.specific_path,'*.txt')
        #                 specific_file_list = glob.glob(specific_file_path)
        specific_file_list = glob.glob(sam_path)

        total_file = len(specific_file_list)
        print(f'{total_file} target file list have been created!')
        print('start concate file...')
        t = 0
        for n in specific_file_list:
            df = pd.read_csv(self.result_file, index_col=False)
            sample_name = n.split('/')[-1].split('.')[0]
            dfn = pd.read_csv(
                n, sep='\t', skiprows=[
                    0, 1, 2, 3, 4, 5, 6, 7, 8, 9
                ]).rename(columns={'All_Reads.fr_count': sample_name})
            dfm = pd.merge(df, dfn, how='left', on='insert_size').fillna(0)
            dfm.to_csv(self.result_file, index=False)
            t += 1
        if t == total_file:
            print('Concated done,ready to cal median insertszie!')
        inpath = os.path.join(self.sample_path, "fragment_study")
        infile = file_name
        median_isnertsize(inpath, infile)

    def runAll(self):
        #self.create_insertsize_mtr()
        self.file_concate()


if __name__ == '__main__':
    sample_path = sys.argv[1]
    #target_path = 'fragment_study/'
    #specific_path = 'region_splited_sam_insertsize'
    #result_path = '/mnt/hdd/usr/jiaohuanmin/home/jiaohm/data/liuy/fragment_study/LLI/insertsize'
    file_name = 'mis_10_insertsize_mtr.csv'
    file_concate = FileConcate(sample_path, file_name)
    file_concate.runAll()