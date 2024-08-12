# MEFI_code
# Aberrant fragmentomic features of circulating cell-free mitochondrial DNA as novel biomarkers in cancer patients
# Contents
  - [Overview](#overview)
  - [System Requirements](#system-requirements)
  - [Installation Guide](#installation-guide)
  - [Demo](#demo)
  - [Results](#results)
  - [License](#license)
  - [Issues](#issues)
  - [Citation](#citation)

# Overview

MEFI_code is a package containing tools for comprehensive analysis of ccf-mtDNA fragmentomic features, that is capable of dealing with sequencing data. MEFI_code also can be used to calculate the probability of developing cancer and to trace its origin.

# System Requirements

**Hardware Requirements**

MEFI_code requires only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB
CPU: 4+ cores, 3.3+ GHz/core

**Software Requirements**

*OS Requirements*

This package is supported for MacOS,Windos and Linux systems.The package development version is tested on Linux operating systems on CentOS 7.9.

The CRAN package should be compatible with Windows, Mac, and Linux operating systems.

Before setting up the MEFI_code package, users should have R version 3.4.0 or higher , randomForest and several packages set up from CRAN should installed.Python scripts is tested on 3.10.9 and python3 is recommended. To avoid some conflicts between softwares, anaconda or miniconda package manager and virtual environment are also recommended.

*Packages and Softs Dependencies*

MEFI_code mainly depends on the Python scientific stack.
```
numpy==1.23.5
pandas==1.5.3
scipy==1.10.0
pysam==0.22.0
samtools==1.15.1
picard-tools==1.81
```

# Installation Guide

1.Install anaconda
2.Create and activate virtul envirment
```
conda create -n your_env_name python=3.10.9
source activate your_env_name
pip install pandas==1.5.3 scipy==1.10.0 pysam==0.22.0
conda install samtools
conda install r-base=4.3.0
#R language mode insatll model package
install.packages ("randomForest")
```

# Demo
Two test sample are under example root and all the features exhibited in the paper could be computed by running the shell script run_fragment_study.sh as follow once:
```
bash run_fragment_study.sh
```

# Results

There are 12 directories and 30 files under the example/fragment_study as the tree structure:
```
fragment_study
├── fragment_study/insertszie
│   ├── fragment_study/insertszie/cumsumfrequency
│   │   ├── fragment_study/insertszie/cumsumfrequency/mis10_insertsize_cufrequency.csv
│   │   ├── fragment_study/insertszie/cumsumfrequency/mis10_insertsize_cumfrequency50%_insertsize.csv
│   │   └── fragment_study/insertszie/cumsumfrequency/mis10_insertsize_excluded_sample.csv
│   └── fragment_study/insertszie/mis_10_insertsize_mtr.csv
├── fragment_study/mis_10_endmotifs
│   ├── fragment_study/mis_10_endmotifs/endmotifs_50_150.csv
│   ├── fragment_study/mis_10_endmotifs/endmotifs_L150.csv
│   ├── fragment_study/mis_10_endmotifs/endmotifs_total.csv
│   └── fragment_study/mis_10_endmotifs/MDS
│       └── fragment_study/mis_10_endmotifs/MDS/four_bases_enmotifs_entropy.csv
└── fragment_study/mis_10_splitedBy125
    ├── fragment_study/mis_10_splitedBy125/125_splited_sam_file
    │   ├── fragment_study/mis_10_splitedBy125/125_splited_sam_file/LB-116cf.mis.10.L125.sam
    │   ├── fragment_study/mis_10_splitedBy125/125_splited_sam_file/LB-116cf.mis.10.S125.sam
    │   ├── fragment_study/mis_10_splitedBy125/125_splited_sam_file/LB-120cf.mis.10.L125.sam
    │   └── fragment_study/mis_10_splitedBy125/125_splited_sam_file/LB-120cf.mis.10.S125.sam
    ├── fragment_study/mis_10_splitedBy125/depth_file
    │   ├── fragment_study/mis_10_splitedBy125/depth_file/LB-116cf.mis.10.L125.sam.txt
    │   ├── fragment_study/mis_10_splitedBy125/depth_file/LB-116cf.mis.10.S125.sam.txt
    │   ├── fragment_study/mis_10_splitedBy125/depth_file/LB-120cf.mis.10.L125.sam.txt
    │   └── fragment_study/mis_10_splitedBy125/depth_file/LB-120cf.mis.10.S125.sam.txt
    ├── fragment_study/mis_10_splitedBy125/depth_stat_LS125_cunt.csv
    ├── fragment_study/mis_10_splitedBy125/depth_stat_LS125_ratio.csv
    ├── fragment_study/mis_10_splitedBy125/depth_stat_LS125_zfsd.csv
    └── fragment_study/mis_10_splitedBy125/FSD
        ├── fragment_study/mis_10_splitedBy125/FSD/AUDC
        │   ├── fragment_study/mis_10_splitedBy125/FSD/AUDC/255_audc_nega.csv
        │   └── fragment_study/mis_10_splitedBy125/FSD/AUDC/255_audc_posi.csv
        ├── fragment_study/mis_10_splitedBy125/FSD/corr
        │   └── fragment_study/mis_10_splitedBy125/FSD/corr/after_smothed_corr_with_mean_of_baseline.csv
        ├── fragment_study/mis_10_splitedBy125/FSD/ED
        │   ├── fragment_study/mis_10_splitedBy125/FSD/ED/euclidean_of_zfsd_and_HC_median.csv
        │   ├── fragment_study/mis_10_splitedBy125/FSD/ED/euclidean_of_zfsd_and_HC_median_in_255.csv
        │   └── fragment_study/mis_10_splitedBy125/FSD/ED/excluded_samples.csv
        └── fragment_study/mis_10_splitedBy125/FSD/peaks
            ├── fragment_study/mis_10_splitedBy125/FSD/peaks/diff2HC_peaks_feature.csv
            ├── fragment_study/mis_10_splitedBy125/FSD/peaks/similar2HC_peaks_feature.csv
            ├── fragment_study/mis_10_splitedBy125/FSD/peaks/similar_diff_peak_counts.csv
            ├── fragment_study/mis_10_splitedBy125/FSD/peaks/smothed_para_51_1_filteredBywidth_5_peaks_feature.csv
            └── fragment_study/mis_10_splitedBy125/FSD/peaks/zfsd_smothed_para_51_1.csv
```


# License

This project is covered under the Apache 2.0 License.

# Issues

If you're having trouble, notice a bug, or want to contribute (such as a fix to the bug you may have just found) feel free to open a git issue or pull request. Enjoy!









