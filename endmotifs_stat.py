import os
import sys
import glob
import pandas as pd
from collections import Counter
from functools import reduce
import pysam
from scipy.stats import entropy
import math

from reverse_complementary_DNA import complementary

def get_endmotifs(sample_path, base_num, target_path, end_total, MDS,MDS_path,
                  end_50_150=None, end_L150=None):
    #sam_files = glob.glob(os.path.join(sample_path, '*/*.mis.10.sam'))
    std_AT = 0.2782  # endmotifs std A and T base ratio in mtDNA reference,mis10AT=0.2782,GC=0.2218
    std_GC = 0.2218  # endmotifs std G and C base ratio in mtDNA reference
    n = int(base_num)
    #endmotifs base type generating
    base_list = ['A', 'T', 'G', 'C']
    motif_set_n = dict((k, []) for k in range(1, n + 1))
    motif_set_n[1] = base_list
    for k in range(2, n + 1):
        motif_set_n[k] = [i + j for i in base_list for j in motif_set_n[k - 1]]

    sample_list = []
    list_total = []
    f = open(os.path.join(sample_path.strip(), 'sample_name.txt'))
    sam_files = [
        os.path.join(sample_path, i.strip(),
                     i.strip() + '.mis.10.bam') for i in f.readlines()
    ]
    f.close()
    exc = []
    for i in sam_files:
        if os.path.getsize(i) > 1:  #check bam file is not empty
            name = os.path.split(i)[1].split('.')[0]
            sample_list.append(name)
            seqDic_t = dict((k, []) for k in range(1, n + 1))
            #insertsize = []
            sam_ali = pysam.AlignmentFile(i, 'rb')  #open sam
            for read in sam_ali.fetch():
                qaseq = read.query_alignment_sequence
                flag = read.flag
                if flag == 83:
                    for k in range(1, n + 1):
                        seqDic_t[k].append(complementary(qaseq[-k:]))
                if flag == 99:
                    for k in range(1, n + 1):
                        seqDic_t[k].append(qaseq[:k])
            sam_ali.close()

            total_count = {
                m: seqDic_t[k].count(m)
                for k in range(1, n + 1) for m in motif_set_n[k]
            }
            df_t = pd.DataFrame(total_count, index=[0])
            list_total.append(df_t)
        else:
            name = os.path.split(i)[1].split('.')[0]
            print(f"{name} mis 10 empty")
    df_total = pd.concat(list_total, axis=0)
    df_total = df_total.reset_index()
    #endmotifs ratio starded by ref base ratio
    ATCG_count = df_total.loc[:, ['A', 'T', 'G', 'C']].sum(axis=1)
    df_total['A_std%'] = df_total['A'] / (ATCG_count * std_AT)
    df_total['T_std%'] = df_total['T'] / (ATCG_count * std_AT)
    df_total['G_std%'] = df_total['G'] / (ATCG_count * std_GC)
    df_total['C_std%'] = df_total['C'] / (ATCG_count * std_GC)
    df_total = df_total.drop('index', axis=1)
    df_total.insert(0, 'sample', sample_list)

    #results saving
    end_outpath = os.path.join(sample_path, target_path)
    os.makedirs(end_outpath, exist_ok=True)
    total_out_path = os.path.join(end_outpath, end_total)
    df_total.to_csv(total_out_path, index=False)
    print(end_total + ' processed done,file saved!')

    #compute MDS
    end = df_total[motif_set_n[n]]
    end = end.apply(lambda x: x / x.sum(), axis=1)
    df_total_r = df_total.iloc[:,1:-4].apply(lambda x:x/x.sum(),axis=0)
    df_total_r = pd.concat([df_total[["sample"]],df_total_r],axis=1)
    df_total_r.to_csv(os.path.join(end_outpath,"endmotifs_total_ratio.csv"), index=False)
    
    pk = end.T
    sub = pd.DataFrame({
        'sample': df_total['sample'],
        'entropy': entropy(pk, base=math.pow(4,n))
    })
    os.makedirs(os.path.join(sample_path, MDS_path), exist_ok=True)
    sub.to_csv(os.path.join(sample_path, MDS_path,str(base_num)+'_'+MDS))
    print("MDS saved!")

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('sample_path', type=str, help='Path to the sample directory')
    parser.add_argument('base_num', type=int, help='Number of bases to consider')
    args = parser.parse_args()

    sample_path = args.sample_path
    base_num = args.base_num
    end_total = 'endmotifs_total.csv'
    MDS = 'bases_enmotifs_entropy.csv'
    MDS_path = "fragment_study/mis_10_endmotifs_4base/MDS"
    
    get_endmotifs(sample_path=sample_path, 
                  base_num=base_num, 
                  target_path='fragment_study/mis_10_endmotifs_4base', 
                  end_total=end_total, 
                  MDS=MDS,
                  MDS_path=MDS_path)
