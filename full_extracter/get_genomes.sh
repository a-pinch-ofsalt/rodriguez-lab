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
