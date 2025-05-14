#!/bin/bash
#BSUB -P acc_oscarlr
#BSUB -q premium
#BSUB -n 4
#BSUB -W 00:05
#BSUB -o pggb_yes.out


# git clone https://github.com/ekg/guix-genomics
# cd guix-genomics
# echo HOWDY
# pwd
# GUIX_PACKAGE_PATH=. guix package -i pggb




module load anaconda3 samtools
# source 

# function create_env {
#     conda init
#     conda create -n pggb-env -c bioconda -c conda-forge odgi smoothxg wfmash=0.13.1 -y
#     conda activate pggb-env
#     conda install -c bioconda pggb -y
# }

# create_env

# conda activate pggb-env



data_folder="/sc/arion/work/hiciaf01/data/2025-05-12_ref_genomes"
small_data_folder="/sc/arion/work/hiciaf01/results/2025-05-12_extracting_igh"

cat ${small_data_folder}/hg19_region.fa ${small_data_folder}/hg38_region.fa > combined.fa
samtools faidx combined.fa

echo OOGABOOGA
pggb -i combined.fa -o output -s 5000 -p 95 -n 2 -t 1
