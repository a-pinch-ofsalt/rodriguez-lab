#!/bin/bash
#BSUB -J pggb
#BSUB -n 4
#BSUB -R "rusage[mem=8000]"
#BSUB -W 4:00
#BSUB -o pggb.out
#BSUB -e pggb.err
#BSUB -P testing_pggb

git clone https://github.com/ekg/guix-genomics
cd guix-genomics
GUIX_PACKAGE_PATH=. guix package -i pggb


data_folder="/sc/arion/work/hiciaf01/data/2025-05-12_ref_genomes"
small_data_folder="/sc/arion/work/hiciaf01/results/2025-05-12_extracting_igh"

cat ${small_data_folder}/hg19_region.fa ${small_data_folder}/hg38_region.fa > combined.fa
samtools faidx combined.fa

pggb -i combined.fa -o output -s 5000 -p 95 -n 2 -t 1
