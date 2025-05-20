#!/bin/bash

hg38_file_general_match=$(ls *hg38.fa* 2>/dev/null | head -n 1)
hg19_file_general_match=$(ls *hg19.fa* 2>/dev/null | head -n 1)

if [ ! -n "${hg38_file_general_match}" ]; then
	echo "Downloading hg38.fa.gz"
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
fi

if [ ! -n "${hg19_file_general_match}" ]; then
	echo "Downloading hg19.fa.gz"
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
fi

if [ ! -f hg38.fa ]; then
	echo "Unzipping hg38.fa.gz"
	gunzip hg38.fa.gz	
fi

if [ ! -f hg19.fa ]; then
	echo "Unzipping hg19.fa.gz"
	gunzip hg19.fa.gz
fi

function minigraph {
    
    git clone https://github.com/lh3/minigraph

    cd minigraph && make

    ./minigraph -cxggs -t16 ../hg19_region.fa ../hg38_region.fa > minigraph_test.gfa
}
function run_pggb {
    
    module load anaconda3 samtools


    function create_env {
	conda init
	conda create -n pggb-env -c bioconda -c conda-forge odgi smoothxg wfmash=0.13.1 -y
	conda activate pggb-env
	conda install -c bioconda pggb -y
    }
    
    create_env
    
    conda activate pggb-env
    
    cat hg19_region.fa hg38_region.fa > combined.fa
    samtools faidx combined.fa
    
    pggb -i combined.fa -o output -s 5000 -p 95 -n 2 -t 1
}
function cut_seqs {
    if ! command -v samtools &> /dev/null; then
	echo "Installing samtools..."
	module load samtools
    else
	echo "Samtools already installed..."
    fi
    
    hg38_region="chr14:105039756-107043718"
    hg19_region="chr14:105505611-107349540"
    
    samtools faidx hg38.fa ${hg38_region} > hg38_region.fa
    samtools faidx hg19.fa ${hg19_region} > hg19_region.fa
}

# download hg_gencode_v47.bed and hg19_vega.bed from https://genome.ucsc.edu/cgi-bin/hgTables
function backup_genes {

    cp hg38_gencode_v47.bed hg38_gencode_v47_backup.bed
    cp hg19_vega.bed hg19_vega_backup.bed
}
function filter_ighv_genes {
    grep 'IGHV' hg38_gencode_v47_backup.bed | grep -v 'alt' > temp_hg38.bed
    mv temp_hg38.bed hg38_gencode_v47_backup.bed

    grep 'IGHV' hg19_vega_backup.bed | grep -v 'alt' > temp_hg19.bed
    mv temp_hg19.bed hg19_vega_backup.bed
}
function subtract_start_position_and_set_proper_labels {
    awk -v start_pos=105039755 '{ $1 = "chr14:105039756-107043718"; $2 -= start_pos; $3 -= start_pos; print }' hg38_gencode_v47_backup.bed > temp_hg38.bed
    mv temp_hg38.bed hg38_gencode_v47_backup.bed
    awk -v start_pos=105505611 '{ $1 = "chr14:105505611-107349540"; $2 -= start_pos; $3 -= start_pos; print }' hg19_vega_backup.bed > temp_hg19.bed
iead hg19_vega.bed
    echo "TEST"
    head temp_hg19.bed
    mv temp_hg19.bed hg19_vega_backup.bed
}
function inject_genes {
    odgi inject -i pggb.og -b hg38_gencode_v47_backup.bed -o pggb_hg38_ighv.og
    odgi inject -i pggb_hg38_ighv.og -b hg19_vega_backup.bed -o pggb_injected.og
}

function get_genes_in_SVs {
    odgi position -i pggb.og -p IGHV3-64D
}

#backup_genes
#filter_ighv_genes
#subtract_start_position_and_set_proper_labels
inject_genes

