#!/bin/bash

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
