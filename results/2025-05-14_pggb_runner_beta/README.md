
# Most recent mini-projct: PGGB on Minerva
PGGB is a difficult to setup on Minerva. Here's how you can do it:

## How to run PGGB on Minerva using this repo
1. `cd` into `results/2025-05-14_pggb_runner_beta`.
2. Run `get_genomes.sh`, `slice_and_index_genomes.sh`, and `pggb_runner.sh`, in that order.
3. In the `output`directory, you should see files ending in `.gfa`, `.pfa`, `.log`, and `.yml`. If `.gfa` is missing, check the `.log`.

##  How to setup PGGB on Minerva yourself

Read the documentation [https://pggb.readthedocs.io/en/latest/]
1. Use bioconda
2. Use whmash 0.13.1

## How to view output pangenome with Bandage 
1. Download Bandage [https://rrwick.github.io/Bandage/]
2. Click `File > Load graph`
3. Locate the `.gfa`
4. Click `Draw graph`

## Known possible errors
1. Problem: `Argument 'INT' received invalid value type '0.001'`. Solution: try version 0.13.1 of whmash `conda install -c bioconda wfmash=0.13.1`.
