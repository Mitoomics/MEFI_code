import os, sys
from collections import Counter
from functools import reduce
import numpy as np
import pandas as pd
from scipy.signal import chirp, find_peaks, peak_widths, savgol_filter
from scipy.spatial import distance

from window_cut import cut

# 1 good data 个人参考基因组比对，提取mtDNA reads
# 2 质控：mis10 过滤
# 3 以125为标准划分mis10为长片段（L125）和短片段（S125）
# 4 计算长片段和短片段depth
# 5 fsd = Ldepth/Sdepth
# 6 fsd数据清洗（缺失值填充，inf值填充，sklearn或者pandas方法填充）
# 7 fsd数据标准化（z-score）
# 8 zfsd 曲线特征计算
## 曲线峰特征（曲线平滑SciPy）
## 曲线下面积
## 曲线与baseline 欧几里得距离（SciPy）
## 曲线与HC曲线相关性（raw zfsd data,平滑后相关性）


# 计算实验样本平滑后曲线峰特征
class ZFSD_Feature_Cal:
    def __init__(
        self,
        root_dir,
        zfsd_path = "fragment_study/mis_10_splitedBy125/depth_stat_LS125_zfsd.csv",
        baseline_zfsd="/mnt/hdd/usr/jiaohuanmin/home/jiaohm/data/liuy/fragment_study/control_data/zfsd/HC_sample_FSD_with_mean_median.csv",
        baseline_smothed="/mnt/hdd/usr/jiaohuanmin/home/jiaohm/data/liuy/fragment_study/LLI/FSD/control/30_HC_sample_zfsd_smothed.csv",
        baseline_peaks_feature="/mnt/hdd/usr/jiaohuanmin/home/jiaohm/data/liuy/fragment_study/LLI/original_data/30_HC_mean_data_after_smothed_para_51_1_peaks_filteredBywidth_5_peaks_feature.csv"
    ):
        self.root_dir = root_dir
        self.zfsd_path = os.path.join(root_dir, zfsd_path)
        self.baseline_zfsd = baseline_zfsd
        self.baseline_smothed = baseline_smothed
        self.baseline_peaks_feature = baseline_peaks_feature
        os.chdir(self.root_dir)
        print(os.getcwd())

    def peaks(self,
              peaks_path="fragment_study/mis_10_splitedBy125/FSD/peaks",
              corr_path="fragment_study/mis_10_splitedBy125/FSD/corr",
              out1="smothed_para_51_1_filteredBywidth_5_peaks_feature.csv",
              out2="zfsd_smothed_para_51_1.csv",
              out3="similar2HC_peaks_feature.csv",
              out4="diff2HC_peaks_feature.csv",
              out5="similar_diff_peak_counts.csv",
              out6="after_smothed_corr_with_mean_of_baseline.csv"):
        baseline_peaks = pd.read_csv(self.baseline_peaks_feature,
                                     index_col=False)
        smothed = []
        feature_exp = []
        feature_s = []
        feature_n = []
        peaks_cunt = {}
        s_cunts = {}
        d_cunts = {}
        dft = pd.read_csv(self.zfsd_path, index_col=False)
        for s in dft.columns:
            name = s.split("_")[1].split(".")[0]
            y = savgol_filter(dft[s], 51, 1,
                              mode="nearest")  # 根据均值确定平滑曲线参数第一个参数51
            peaks_positive, _ = find_peaks(y, width=5)
            peaks_negative, _ = find_peaks(-y, width=5)
            peaks_cunt[name] = len(peaks_positive) + len(peaks_negative)
            w_positive = peak_widths(y, peaks_positive, rel_height=1)
            w_negative = peak_widths(-y, peaks_negative, rel_height=1)
            sub_f = pd.DataFrame({name: y})
            smothed.append(sub_f)
            sub_chara_p = pd.DataFrame({
                "pos": peaks_positive,
                name + "_peakvalue": y[peaks_positive],
                name + "_peakwidth": w_positive[0],
                name + "_widthhight": w_positive[1],
                name + "_widthleft": w_positive[2],
                name + "_widthright": w_positive[3],
            })
            sub_chara_n = pd.DataFrame({
                "pos": peaks_negative,
                name + "_peakvalue": y[peaks_negative],
                name + "_peakwidth": w_negative[0],
                name + "_widthhight": w_negative[1],
                name + "_widthleft": w_negative[2],
                name + "_widthright": w_negative[3],
            })
            sub_feature_exp = (pd.concat(
                [sub_chara_p, sub_chara_n],
                axis=0).sort_values(by="pos").reset_index())
            sub_feature_exp = sub_feature_exp.drop("index", axis=1)
            feature_exp.append(sub_feature_exp)
            similar = set()
            for i in sub_feature_exp["pos"]:
                for j in range(21):
                    if i + j in list(baseline_peaks["pos"]):
                        similar.add(i)
                    if abs(i - j) in list(baseline_peaks["pos"]):
                        similar.add(i)

            dfs = sub_feature_exp.loc[sub_feature_exp["pos"].isin(similar)]
            dfs = dfs.rename(columns={"pos": "spos"})
            feature_s.append(dfs)
            s_cunts[name] = len(similar)

            var = set(sub_feature_exp["pos"]) - similar
            dfn = sub_feature_exp.loc[sub_feature_exp["pos"].isin(var)]
            dfn = dfn.rename(columns={"pos": "dpos"})
            feature_n.append(dfn)
            d_cunts[name] = len(var)

        data_peaks = reduce(lambda i, j: pd.merge(i, j, on="pos", how="outer"),
                            feature_exp)
        os.makedirs(peaks_path, exist_ok=True)
        data_peaks.to_csv(os.path.join(self.root_dir, peaks_path, out1),
                          index=False)
        data_filtered = pd.concat(smothed, axis=1)
        data_filtered.to_csv(os.path.join(self.root_dir, peaks_path, out2),
                             index=False)

        data_s = reduce(lambda i, j: pd.merge(i, j, on="spos", how="outer"),
                        feature_s)
        data_n = reduce(lambda i, j: pd.merge(i, j, on="dpos", how="outer"),
                        feature_n)
        data_s.to_csv(os.path.join(self.root_dir, peaks_path, out3),
                      index=False)
        data_n.to_csv(os.path.join(self.root_dir, peaks_path, out4),
                      index=False)

        df_s = (pd.DataFrame.from_dict(
            s_cunts, orient="index",
            columns=["similar_peak_cunts"
                     ]).reset_index().rename(columns={"index": "sample"}))
        df_n = (pd.DataFrame.from_dict(
            d_cunts, orient="index",
            columns=["different_peak_cunts"
                     ]).reset_index().rename(columns={"index": "sample"}))
        df_cunts = (pd.DataFrame.from_dict(
            peaks_cunt, orient="index",
            columns=["peak_counts"
                     ]).reset_index().rename(columns={"index": "sample"}))
        cunt_list = [df_s, df_n, df_cunts]
        df_sn_cunt = reduce(lambda x, y: pd.merge(x, y, on="sample"),
                            cunt_list)
        df_sn_cunt.to_csv(os.path.join(self.root_dir, peaks_path, out5),
                          index=False)

        # 计算实验样本与baseline样本平滑后zfsd曲线相关性
        base_smothed = pd.read_csv(self.baseline_smothed,
                                   usecols=["Mean"],
                                   index_col=False)
        os.makedirs(corr_path, exist_ok=True)
        pearson = {}
        for col in data_filtered.columns:
            pearson[col] = base_smothed["Mean"].corr(data_filtered[col])
        pd.DataFrame.from_dict(pearson, orient="index",
                               columns=["pearson"]).to_csv(
                                   os.path.join(self.root_dir, corr_path,
                                                out6))
        print("peaks feature stat done!")
        print("person corr results saved!")

    def euclidean_distance(
            self,
            ed_path="fragment_study/mis_10_splitedBy125/FSD/ED",
            output1="euclidean_of_zfsd_and_HC_median.csv",
            output2="euclidean_of_zfsd_and_HC_median_in_255.csv",
            output3="excluded_samples.csv"):
        # zfsd曲线与参考df_h相似度：欧几里得距离
        df_h = pd.read_csv(self.baseline_zfsd, index_col=False)
        dft = pd.read_csv(self.zfsd_path, index_col=False)
        df_255_list = []
        df_total_list = []
        sample_excl = []
        ou_dic_total = {}
        ou_dic_255 = dict((k, []) for k in dft.columns)
        for col in dft.columns:
            try:
                if not dft[col].empty:
                    # 计算整条曲线与baseline相似度
                    ou_dic_total[col] = distance.euclidean(
                        dft[col], df_h["median"])
                    # 计算255个区间内曲线相似度（欧氏距离）
                    base_cut_list = cut(df_h["median"], 65)
                    sample_cut_list = cut(dft[col], 65)
                    for j in range(255):
                        oui = distance.euclidean(sample_cut_list[j],
                                                 base_cut_list[j])
                        ou_dic_255[col].append(oui)
                else:
                    sample_excl.append(col)
                    del ou_dic_255[col]
                    pass
            except:
                print(col)
                sample_excl.append(col)
                ou_dic_255.pop(col)
        df_ou_255 = pd.DataFrame(ou_dic_255)
        # df_255_list.append(df_ou_255)
        df_ou_total = pd.DataFrame({
            "sample": ou_dic_total.keys(),
            "Euclidean": ou_dic_total.values()
        })
        # df_55 = pd.concat(df_255_list, axis=1)
        os.makedirs(ed_path, exist_ok=True)
        df_ou_total.to_csv(os.path.join(self.root_dir, ed_path, output1),
                           index=False)
        df_ou_255.to_csv(os.path.join(self.root_dir, ed_path, output2),
                         index=False)
        df_ex = pd.DataFrame({"sample_excluded": sample_excl})
        df_ex.to_csv(os.path.join(self.root_dir, ed_path, output3),
                     index=False)
        print("euclidean_distance calculate done,results saved~")

    def audc(self,
             audc_path="fragment_study/mis_10_splitedBy125/FSD/AUDC",
             output1="255_audc_posi.csv",
             output2="255_audc_nega.csv"):
        # 计算曲线下面积
        dft = pd.read_csv(self.zfsd_path, index_col=False)
        audc_dic_255_pos = dict((k, []) for k in dft.columns)
        audc_dic_255_nega = dict((k, []) for k in dft.columns)
        for s in dft.columns:
            sample_cut_list = cut(dft[s], 65)
            for j in range(255):
                cuti = np.array(sample_cut_list[j])
                audci_pos = cuti[cuti >= 0].sum()
                audci_nega = cuti[cuti < 0].sum()
                audc_dic_255_pos[s].append(audci_pos)
                audc_dic_255_nega[s].append(audci_nega)

        # results stat
        df_audc_dic_255_pos = pd.DataFrame(audc_dic_255_pos)
        df_audc_dic_255_nega = pd.DataFrame(audc_dic_255_nega)
        os.makedirs(audc_path, exist_ok=True)
        df_audc_dic_255_pos.to_csv(os.path.join(self.root_dir, audc_path,
                                                output1),
                                   index=False)
        df_audc_dic_255_nega.to_csv(os.path.join(self.root_dir, audc_path,
                                                 output2),
                                    index=False)
        print("audc calculated done! results saved")

    def run(self):
        print("compute peaks!")
        self.peaks()
        print("compute euclidean distance")
        self.euclidean_distance()
        print("cal audc!")
        self.audc()


if __name__ == "__main__":
    root_dir = sys.argv[1]
    #root_dir = "/mnt/data1/Seq_data/20230901/HCS_Plasma_DOU_mtDNA-CAP/Analysis_private-mtDNA_V1"
    fsd_feature_get = ZFSD_Feature_Cal(root_dir)
    fsd_feature_get.run()