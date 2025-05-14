git clone https://github.com/lh3/minigraph

data_folder="/sc/arion/work/hiciaf01/results/2025-05-12_extracting_igh"

cd minigraph && make

./minigraph -cxggs -t16 ${data_folder}/hg19_region.fa ${data_folder}/hg38_region.fa > lol.gfa
