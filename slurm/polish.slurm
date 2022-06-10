#!/bin/bash
#!

#SBATCH -J POLISH
#SBATCH -A <<REPLACE>>
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=36:00:00
#SBATCH --mail-type=FAIL
#SBATCH --exclusive
#SBATCH -p <<REPLACE>>

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')

#! Load whichever modules you need
. /etc/profile.d/modules.sh                # This line enables the module command
module purge                               # Removes all modules still loaded

if [ $build = "GRCh37" ]; then
 buildInfix="hg19"
 ref=$build_37
else
 buildInfix="hg38"
 ref=$build_38
fi

aligner="minimap"

#! Source config file
source $1

#! Enable conda usage
source $conda_path

conda activate pysam

for caller in "cutesv" "svim" ; do

    echo "Processing $sample_name $caller"

    bamFile="${out_dir}/alignments/${build}/${aligner}/${sample_name}.${buildInfix}.bam"

    python3 /path/to/getGoodAlts.py \
    ${out_dir}/variant_calling/GRCh38/${caller}/${aligner}/$sample_name.${buildInfix}.${aligner}.${caller}.min3.vcf.gz \
    $bamFile \
    /path/to/polish/$sample_name.${buildInfix}.${aligner}.${caller}.min3.polish.vcf

done

conda deactivate
