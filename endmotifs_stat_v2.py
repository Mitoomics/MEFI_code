import pandas as pd
import os, sys
import glob
from collections import Counter
from functools import reduce
import pysam
from reverse_complementary_DNA import complementary
from scipy.stats import entropy


def get_endmotifs(sample_path, base_num, target_path, end_total, MDS,
                  end_50_150, end_L150):
    #sam_files = glob.glob(os.path.join(sample_path, '*/*.mis.10.sam'))
    std_AT = 0.2782  # endmotifs std A and T base ratio in mtDNA reference,mis10AT=0.2782,GC=0.2218
    std_GC = 0.2218  # endmotifs std G and C base ratio in mtDNA reference
    n = base_num
    #endmotifs base type generating
    base_list = ['A', 'T', 'G', 'C']
    motif_set_n = dict((k, []) for k in range(1, n + 1))
    motif_set_n[1] = base_list
    for k in range(2, n + 1):
        motif_set_n[k] = [i + j for i in base_list for j in motif_set_n[k - 1]]

    sample_list = []
    list_total = []
    list_50_150 = []
    list_L150 = []
    #sam_files = glob.glob(os.path.join(sample_path,"fragment_study/region_splited_sam_7sandother","*.mis.10.other.sam"))
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
            insertsize = []
            sam_ali = pysam.AlignmentFile(i, 'rb')  #open sam
            for read in sam_ali.fetch():
                qaseq = read.query_alignment_sequence
                flag = read.flag
                if flag == 83:
                    insert = read.template_length
                    insertsize.append(abs(insert))
                    for k in range(1, n + 1):
                        seqDic_t[k].append(complementary(qaseq[-k:]))
                if flag == 99:
                    insert = read.template_length
                    insertsize.append(abs(insert))
                    for k in range(1, n + 1):
                        seqDic_t[k].append(qaseq[:k])
            sam_ali.close()

            #counts of reads in every base type
            total_count = {
                m: seqDic_t[k].count(m)
                for k in range(1, n + 1) for m in motif_set_n[k]
            }

            seqDic_sub1 = dict((k, []) for k in range(1, n + 1))
            seqDic_sub2 = dict((k, []) for k in range(1, n + 1))
            for s in range(len(insertsize)):
                if (insertsize[s] >= 50 and insertsize[s] <= 150):
                    for k in range(1, n + 1):
                        seqDic_sub1[k].append(seqDic_t[k][s])
                if (insertsize[s] > 150):
                    for k in range(1, n + 1):
                        seqDic_sub2[k].append(seqDic_t[k][s])

            sub1_count = {
                m: seqDic_sub1[k].count(m)
                for k in range(1, n + 1) for m in motif_set_n[k]
            }
            sub2_count = {
                m: seqDic_sub2[k].count(m)
                for k in range(1, n + 1) for m in motif_set_n[k]
            }
            df_t = pd.DataFrame(total_count, index=[0])
            list_total.append(df_t)
            df_sub1 = pd.DataFrame(sub1_count, index=[0])
            list_50_150.append(df_sub1)
            df_sub2 = pd.DataFrame(sub2_count, index=[0])
            list_L150.append(df_sub2)
        else:
            name = os.path.split(i)[1].split('.')[0]
            print(f"{name} mis 10 empty")
            exc.append(name)
    #with open('/data/usr/jiaohm/data/LLI/excluded_sample/excluded_new_sample.txt',mode='a') as f:
    #    f.writelines([i+'\n' for i in exc])
    #endmotifs results concat
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
    df_50_150 = pd.concat(list_50_150, axis=0)
    df_50_150 = df_50_150.reset_index()
    df_50_150 = df_50_150.drop('index', axis=1)
    df_50_150.insert(0, 'sample', sample_list)
    df_L150 = pd.concat(list_L150, axis=0)
    df_L150 = df_L150.reset_index()
    df_L150 = df_L150.drop('index', axis=1)
    df_L150.insert(0, 'sample', sample_list)

    #results saving
    end_outpath = os.path.join(sample_path, target_path)
    os.makedirs(end_outpath, exist_ok=True)
    total_out_path = os.path.join(end_outpath, end_total)
    out_50_150_path = os.path.join(end_outpath, end_50_150)
    out_L150_path = os.path.join(end_outpath, end_L150)

    df_total.to_csv(total_out_path, index=False)
    df_50_150.to_csv(out_50_150_path, index=False)
    df_L150.to_csv(out_L150_path, index=False)
    print(end_total + ' processed done,file saved!')
    print(end_50_150 + ' processed done,file saved!')
    print(end_L150 + ' processed done,file saved!')

    #计算MDS
    end = df_total[motif_set_n[4]]
    end = end.apply(lambda x: x / x.sum(), axis=1)
    pk = end.T
    sub = pd.DataFrame({
        'sample': df_total['sample'],
        'entropy': entropy(pk, base=256)
    })
    os.makedirs(os.path.join(sample_path, MDS_path), exist_ok=True)
    sub.to_csv(os.path.join(sample_path, MDS_path, MDS))
    print("MDS saved!")


if __name__ == '__main__':
    sample_path = sys.argv[1]
#     sample_path = '/mnt/data1/Seq_data/20231124/LC_Plasma_DOU_mtDNA-CAP/Analysis_private-DNA_V1'
    base_num = 4
    end_total = 'endmotifs_total.csv'
    MDS = "four_bases_enmotifs_entropy.csv"
    end_50_150 = 'endmotifs_50_150.csv'
    end_L150 = 'endmotifs_L150.csv'
    target_path = 'fragment_study/mis_10_endmotifs'
    MDS_path = "fragment_study/mis_10_endmotifs/MDS"
    get_endmotifs(sample_path, base_num, target_path, end_total, MDS,
                  end_50_150, end_L150)
