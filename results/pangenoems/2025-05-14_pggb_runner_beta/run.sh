#!/bin/bash

function get_ighv_locus {
    if [ ! -n "$(ls *hg38.fa* 2>/dev/null | head -n 1)" ]; then
        wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    fi
    if [ ! -n "$(ls *hg19.fa* 2>/dev/null | head -n 1)" ]; then
        wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
    fi
    gunzip hg38.fa.gz
    gunzip hg19.fa.gz
    module load samtools
    samtools faidx hg38.fa chr14:105039756-107043718 > hg38_region.fa
    samtools faidx hg19.fa chr14:105505611-107349540 > hg19_region.fa
}

function create_pangenome_env {
    #module load anaconda3
    conda init
    conda create -n pangenome-env -c bioconda -c conda-forge pggb odgi smoothxg wfmash=0.13.1 -y

    # function minigraph_create_pangenome_graph {
    #     git clone https://github.com/lh3/minigraph
    #     cd minigraph && make
    # }
}
function create_pangenome_graph {
    conda init
    conda activate pangenome-env
    module load samtools
    cat hg19_region.fa hg38_region.fa > combined.fa
    samtools faidx combined.fa
    pggb -i combined.fa -o output -s 5000 -p 95 -n 2 -t 1

    # ./minigraph -cxggs -t16 ../hg19_region.fa ../hg38_region.fa > minigraph_test.gfa
}
# download hg_gencode_v47.bed and hg19_vega.bed from https://genome.ucsc.edu/cgi-bin/hgTables
function backup_ighv {
    cp hg38_gencode_v47.bed hg38_gencode_v47_backup.bed
    cp hg19_vega.bed hg19_vega_backup.bed
}
function filter_ighv_genes {
    grep 'IGHV' hg38_gencode_v47_backup.bed | grep -v 'alt' > temp_hg38.bed
    mv temp_hg38.bed hg38_gencode_v47_backup.bed

    grep 'IGHV' hg19_vega_backup.bed | grep -v 'alt' > temp_hg19.bed
    mv temp_hg19.bed hg19_vega_backup.bed
}
function subtract_start_position_and_set_proper_labels_for_odgi {
    awk -v start_pos=105039755 '{ $1 = "chr14:105039756-107043718"; $2 -= start_pos; $3 -= start_pos; print }' hg38_gencode_v47_backup.bed > temp_hg38.bed
    mv temp_hg38.bed hg38_gencode_v47_backup.bed
    awk -v start_pos=105505611 '{ $1 = "chr14:105505611-107349540"; $2 -= start_pos; $3 -= start_pos; print }' hg19_vega_backup.bed > temp_hg19.bed
iead hg19_vega.bed
    echo "TEST"
    head temp_hg19.bed
    mv temp_hg19.bed hg19_vega_backup.bed
}

function fix_files {
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' hg38_gencode_v47_backup.bed > hg38_gencode_v47_backup_fixed.bed
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' hg19_vega_backup.bed > hg19_vega_backup_fixed.bed
    mv hg38_gencode_v47_backup_fixed.bed hg38_gencode_v47_backup.bed
    mv hg19_vega_backup_fixed.bed hg19_vega_backup.bed

}

#make sure to run within a job
function inject_genes {
    #bsub -P acc_oscarlr -q interactive -n 8 -W 1 -R span[hosts=1]  -Is /bin/bash
    conda init
    conda activate pangenome-env
    odgi inject -i output/combined.fa.36fc4b3.417fcdf.seqwish.gfa -b hg38_gencode_v47_backup.bed -o pggb_hg38_ighv.og
    #odgi inject -i pggb_hg38_ighv.og -b hg19_vega_backup.bed -o pggb_injected.og
}

function get_genes_in_SVs {
    odgi position -i pggb.og -p IGHV3-64D
}
#download_genomes
#create_pangenome_env
# create_pangenome_graph
# backup_ighv
# filter_ighv_genes
# subtract_start_position_and_set_proper_labels_for_odgi
inject_genes

# fix_files