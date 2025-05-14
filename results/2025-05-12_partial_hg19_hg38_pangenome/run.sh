#!/bin/bash

module load anaconda3 samtools
source ~/.bashrc

function create_env {
    conda init
    conda create -n pggb-env -c bioconda -c conda-forge odgi smoothxg wfmash -y
    conda activate pggb-env
    conda install -c bioconda pggb -y
}
#create_env

conda activate pggb-env

git clone --recursive https://github.com/pangenome/pggb.git

data_folder="/sc/arion/work/hiciaf01/data/2025-05-12_ref_genomes"

cat ${data_folder}/hg19.fa ${data_folder}/hg38.fa > combined.fa
samtools faidx combined.fa

cd pggb
pggb -i ../combined.fa -o output -s 5000 -p 95 -n 2 -t 8
