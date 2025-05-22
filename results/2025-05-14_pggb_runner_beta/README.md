
# Most recent mini-projct: PGGB on Minerva
PGGB is a difficult to setup on Minerva. Here's how you can do it:

use whmash 0.13.1 for pggb 

## Known possible errors
1. Problem: `Argument 'INT' received invalid value type '0.001'`. Solution: try version 0.13.1 of whmash `conda install -c bioconda wfmash=0.13.1`.
2. Using [this env](https://pastebin.com/aauVRSus) I get a corrupted image and misplaced genes,  using [this one](https://pastebin.com/F2qNQmhP) it's fine
![hg38_ighv](https://github.com/user-attachments/assets/65eb0d29-0a78-44c2-85bc-40cd8bfc5b77)

4. You'll have to do some stuff, like pggb pangenome building and odgi injecting in a job, not in the login node.
5. odgi requires this .bed format `chr14:105505611-107349540	1186061	1186592	IGHV3-21`, **tab, not space separated** no header. I tried concatenating two bed files to align ighv genes from two separate genomes, hg38 and hg19, in the same graph, but that's an error.

## Continuing
NA12878: https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Datasets
CHM1.1: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000306695.2/
1. Setup AWS CLI using Conda, use aws s3 ls on the "s3://platinum-pedigree-data/assemblies/NA12878/" folder and cp whatever file you want within there
2. Download CHM1.1 and upload it manually
3. Here I didn't know the start and ends of ighv locus in these genomes, I was thinking BLAST using hg38_region.fa. 
