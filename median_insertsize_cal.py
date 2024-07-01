import os, sys
import pandas as pd


def median_isnertsize(
        inpath,
        infile,
        outfile1="mis10_insertsize_cufrequency.csv",
        outfile2="mis10_insertsize_cumfrequency50%_insertsize.csv",
        excluded="mis10_insertsize_excluded_sample.csv"):
    df = pd.read_csv(os.path.join(inpath,infile), index_col=False)
    merge_list = []
    d = []
    for col in df.iloc[:, 1:].columns:
        a = df[col] / df[col].sum()
        b = a.cumsum()
        df_cumsumfreq = pd.DataFrame({col: b})
        merge_list.append(df_cumsumfreq)
        d.append(a)
    df_cf = pd.concat(merge_list, axis=1)
    df_d = pd.concat(d, axis=1)
    outpath = os.path.join(inpath, "cumsumfrequency")
    os.makedirs(outpath, exist_ok=True)
    df_cf.to_csv(os.path.join(outpath, outfile1), index=False)
    df_d.to_csv(os.path.join(outpath, "denisty.csv"), index=False)

    col_dic = {}
    col_aa = {}
    excluded_sample = []
    for col in df_cf.columns:
        try:
            col_dic[col] = df_cf[col][(df_cf[col] >= 0.5)
                                      & (df_cf[col] < 0.6)].index.to_list()[0]
            col_aa[col] = list(
                df_cf[col][(df_cf[col] >= 0.5) & (df_cf[col] < 0.6)])[0]
        except:
            print(col)
            excluded_sample.append(col)
    df_aa = pd.DataFrame({
        'sample': col_dic.keys(),
        'insertsize': col_dic.values()
    })
    df_bb = pd.DataFrame({'sample': col_aa.keys(), '50%': col_aa.values()})
    dff = pd.merge(df_aa, df_bb, on='sample')
    dff.to_csv(os.path.join(outpath, outfile2), index=False)
    dfe = pd.DataFrame(excluded_sample, columns=['sample'])
    dfe.to_csv(os.path.join(outpath, excluded), index=False)
    print("process done!")


if __name__ == "__main__":
    inpath = sys.argv[1]
    infile = "mis_10_insertsize_mtr.csv"
    median_isnertsize(inpath, infile)
