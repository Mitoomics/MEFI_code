import pandas as pd
import os, sys
import glob
import pysam
from reverse_complementary_DNA import complementary
from functools import reduce


def get_fourebase_end(sample_path, target_path, outpath, endout):
    std_AT = 0.2782  # 末端碱基A、T占比标准化参考值
    std_GC = 0.2218  # 末端C、G碱基标准化参考值
    n = base_num = 4
    #不同长度末端碱基类型统计
    base_list = ['A', 'T', 'G', 'C']
    motif_set_n = dict((k, []) for k in range(1, n + 1))
    motif_set_n[1] = base_list
    for k in range(2, n + 1):
        motif_set_n[k] = [i + j for i in base_list for j in motif_set_n[k - 1]]
    sam_files = glob.glob(
        os.path.join(sample_path, target_path, "*.mis.10.sam"))

    sample_list = []
    list_total = []
    for i in sam_files:
        name = os.path.split(i)[1].rsplit('.', 1)[0]
        print(name)
        sample_list.append(name)
        seqDic_t = dict((k, []) for k in range(1, n + 1))
        insertsize = []
        sam_ali = pysam.AlignmentFile(i, 'r')  #sam文件打开
        for read in sam_ali:
            qaseq = read.query_alignment_sequence
            flag = read.flag
            ref = read.reference_name
            if flag == 83 and ref == 'chrM':
                insert = read.template_length
                for k in range(1, n + 1):
                    seqDic_t[k].append(complementary(qaseq[-k:]))
            if flag == 99 and ref == 'chrM':
                for k in range(1, n + 1):
                    seqDic_t[k].append(qaseq[:k])
        sam_ali.close()

        #统计末端相应末端碱基类型reads数量
        total_count = {
            m: seqDic_t[k].count(m)
            for k in range(1, n + 1) for m in motif_set_n[k]
        }

        df_t = pd.DataFrame(total_count, index=[0])
        list_total.append(df_t)

    #末端结果统计
    dfT = pd.concat(list_total, axis=0)
    dfT = dfT.reset_index()
    df_total = dfT.copy()
    #但末端碱基占比统计并用ref数据标准化
    df_total.loc[:, motif_set_n[4]] = df_total.loc[:, motif_set_n[4]].apply(
        lambda x: x / x.sum(), axis=1)
    ATCG_count = df_total.loc[:, ['A', 'T', 'G', 'C']].sum(axis=1)
    df_total['A_std%'] = df_total['A'] / (ATCG_count * std_AT)
    df_total['T_std%'] = df_total['T'] / (ATCG_count * std_AT)
    df_total['G_std%'] = df_total['G'] / (ATCG_count * std_GC)
    df_total['C_std%'] = df_total['C'] / (ATCG_count * std_GC)
    df_total = df_total.drop('index', axis=1)
    df_total.insert(0, 'sample', sample_list)

    #数据保存
    end_outpath = os.path.join(sample_path, outpath)
    os.makedirs(end_outpath, exist_ok=True)
    total_out_path = os.path.join(end_outpath, endout)
    df_total.to_csv(total_out_path, index=False)
    print(sample_path + ' processed done,file saved!')


if __name__ == "__main__":
    sample_path = sys.argv[1]
    #sample_path = "/data/HCC/02-HC-plasma-and-PBMC/HC-plasma/"
    target_path = "*/"
    outpath = "fragment_study/mis_10_endmotifs"
    endout = 'endmotifs_total.csv'
    get_fourebase_end(sample_path, target_path, outpath, endout)
