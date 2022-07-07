# Introduction

This is a repository meant to accompany the manuscript [TBA]. It contains the scripts that were used to simulate the results, as well as information on the conda environments needed.

To use them, it would be necessary to edit every script and replace the placeholders for paths, which are in the format `/path/to`, and the names for folders, which are `folderName`.

Below, a brief description of each script is provided to make this repository meaningful.

Bash scrips require bash (version `4.2.46(2)-release (x86_64-redhat-linux-gnu)` was used). Python scripts require python >=3.7, and R scripts were executed on R 3.6.

Some of the hevy-load programs were called trhough SLURM. The SLURM scripts include the parameters that were used to call programs that were called that way. For non-SLURM scripts, parameters are included in bash scripts.

Data folder organizations is expected to be:

```
parentFolder
    |- folderName1
    |
    |- folderName2
```

## Abreviations

    - args: command line arguments
    - SURVIVOR: executable to SURVIVOR
    - RE: supporting reads for a structural variant (including insertions)

## Config file

The scripts, particularly the SLURM ones, require a file that can be `source`'d by bash and contains important information, such as paths to bam files, VCF files, etc. These paths can act as output, for aligners and variant callers, or as input, por dowstream analysis. A sample config files is provided.

## Conda environments

A total of four environments are required to reproduce the work in the article (not including R packages).

In the folder `envs`, yml files have been included to recreate these environments. They can be installed in the following manner:

```
conda env create -n env1 --file envs/env1.yml
conda env create -n env2 --file envs/env2.yml
conda env create -n env3 --file envs/env3.yml
conda env create -n env4 --file envs/env4.yml
```
## R dependencies

This is a list of R packages used for data analysis.


# Scripts

## Data fetch and preparation

### prepareNstd152.sh

**Replace**: paths
**Requires**: SURVIVOR, tabix, bgzip, merge.py

This downloads and prepares the VCF file from nstd152 (the HGSVC data).

### prepareNstd215.sh

**Replace**: paths
**Requires**: SURVIVOR, tabix, bgzip, merge.py

This downloads and prepares the VCF file from nstd215 (the Indigen data).  
What a coincidence both files have the same numbers in their name, but in different order.

### prepareOurMe.sh

**Replace**: paths
**Requires**: bcftools

Turn the VCF from the pipeline into a proper VCF.

## Prior stuff

These handle the aligning, variant calling, basic coverage processing, and variants merging.

### slurm/runMinimap2.slurm

**Replace**: SLURM data (as needed)
**Requires**: conda environment env1

Runs the aligner minimap2.

`<<sbatch ...| bash>> runMinimap2.slurm <<config file>>`

### slurm/runNGMLR.slurm

**Replace**: SLURM data (as needed)
**Requires**: conda environment env1

Runs the aligner NGMLR. Only required for reproducing the benchamrk and not necessary for anything else.

`<<sbatch ...| bash>> runNGML.slurm <<config file>>`

### slurm/runLra.slurm

**Replace**: SLURM data (as needed)
**Requires**: conda environment lra

Runs the aligner lra. Only required for reproducing the benchamrk and not necessary for anything else.

`<<sbatch ...| bash>> runLra.slurm <<config file>>`

### slurm/runSniffles.slurm

**Replace**: SLURM data (as needed)
**Requires**: conda environment env1, special version of Sniffles (see next paragraph)

Runs the variant caller Sniffles. It requires version 1.0.12, but since that version had a bug that would generate an invalid VCF header, it needs to be a version of Sniffles 1.0.12 created after September 2020. The bug was fixed on commit "f958698c216a118cd060968e8534e1eeab0ec0bf". Versions available in conda do not include this bugfix (to my knoledge). Sniffles is required for reproducing the benhmark and not necessary for anything else.

`<<sbatch ...| bash>> runSniffles.slurm <<config file>> <<minimap|NGMLR|lra>>`

### slurm/runCuteSV.slurm

**Replace**: SLURM data (as needed)
**Requires**: conda environment env2

Runs the variant caller CuteSV. This caller reports deleted, duplicated, etc. sequences, which can be obtained from the reference and make the VCF take up ~100 GB of space, which is underisable, so those sequences are removed (except for insertions).

`<<sbatch ...| bash>> runCuteSV.slurm <<config file>> <<minimap|NGMLR|lra>>`

### slurm/runSvim.slurm

**Replace**: SLURM data (as needed)
**Requires**: conda environment env2

Runs the variant caller SVIM.

`<<sbatch ...| bash>> runSvim.slurm <<config file>> <<minimap|NGMLR|lra>>`

### slurm/runNanoVar.slurm

**Replace**: SLURM data (as needed)
**Requires**: conda environment env2

Runs the variant caller NanoVar. It also creates a bam file with only reads with MAPQ >= 20, to simulate this option that NanoVar seems to lack, in order to make it simmilar to the other callers. Keep in mind that these program does not seem to be compatible with lra. Only required for reproducing the benchamrk and not necessary for anything else.

`<<sbatch ...| bash>> runNanoVar.slurm <<config file>> <<minimap|NGMLR>>`

### runSurvivor.sh

**Replace**: paths
**Requires**: bcftools

This creates the intra-sample SV consensus.

`bash runSurvivor.sh sampleFolder1 sampleFolder2`

### runSurvivorAll.sh

**Replace**: folderName1, etc. (they can be passed as args); path to SURVIVOR; path to parentFolder
**Requires**: bcftools

This takes the intra-sample merged VCFs (with RE >= 3) and creates inter-sample sets using the lax and strict criteria.

## Benchmark

At this point, it is possible to recreate Supplementary Table 4, with the table benchamrking the combinations of aligner and callers.

### mergeVcfBenchmark

**Replace**: paths
**Requires**: SURVIVOR

This creates merged dataset to facilitate evaluation.

### benchmark.R

**Replace**: paths

This will create an HTML version of Supplementary Table 4, in HTML format.

## TE specific

These scripts are to be executed after the previous ones. They handle the TE analysis.

### getGoodAlts.py

**Requires**: pysam, joblib, Bioconductor, wtdbg2

This, if possible, substitutes the ALT allele produced by the variant callers with a sequence obtained by assembling the surrounding regions.

### slurm/polish.slurm

**Replace**: SLURM parameters, paths
**Requires**: conda, slurm, dependencies for getGoodAlts.py

This was used to call getGoodAlts.py. The python3 command could be used outside SLURM.

### merge.py

**Requires**: pysam, fuzzywuzzy, python-Levhenstein

A script to create a consensus between two or more sets of insertions. Generates an output similar to SURVIVOR. Its key difference with SURVIVOR is that it can take alleles into consideration, and is a bit more flexible with SV length/coordinates.

## mergeIndividual.sh

**Replace**: paths
**Requires**: merge.py

This takes the polished VCFs (with insertions) and creates the intra-sample merge.

## mergeAll.sh

**Replace**: paths
**Requires**: merge.py, RepeatMasker, env3, env4

This takes the intra-sample merged datasets to create a inter-sample merge using merge.py.  
Then, it runs RepeatMasker and creates a bed file with TE insertions.

### slurm/mosdepth.slurm

**Replace**: SLURM parameters, paths
**Requires**: conda env2, slurm

This calculates coverage for TE insertions and deletions. Used for genotyping.

This is run from R.

## getMeDeletions.py

**Requires**: pysam, joblib

This takes a VCF, gets deletions that match a TE in the reference sequence, and creates a new VCF with them.

`python3 getMeDeletions.py /path/to/all.merged2.min3.vcf /path/to/repeatsReferenceTE.bed.gz /paht/to/all.merged2.min3.deletedMeNew.vcf`

## runComparisons.sh

**Replace**: paths
**Requires**: merge.py

This compares the generated dataset against HGSV and Indigen's datasets.

This is run from R.
