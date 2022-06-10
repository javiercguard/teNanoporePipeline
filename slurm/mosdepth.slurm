#!/bin/bash
#!

#SBATCH -J MOSDEPTH
#SBATCH -A <<YOUR_CLUSTER>>
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=36:00:00
#SBATCH --mail-type=FAIL
#SBATCH --exclusive
#SBATCH -p <<YOUR_CLUSTER>>


#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')

#! Load whichever modules you need
. /etc/profile.d/modules.sh                # This line enables the module command
module purge                               # Removes all modules still loaded
module load miniconda/3

## Source the config file
source $1
shift

#! Set ONT enviroment
source $conda_path
mosdepth_env="/path/to/.conda/envs/mosdepth"
conda activate ${mosdepth_env}

#! Define variables
aligner="minimap"

if [ $build = "GRCh37" ]; then
 buildInfix="hg19"
 ref=$build_37
 refBed=${ref_bed_37}
else
 buildInfix="hg38"
 ref=$build_38
 refBed=${ref_bed_38}
fi

bamFile=${out_dir}/alignments/${build}/${aligner}/${sample_name}.${buildInfix}.bam
outputDir=/path/to/retro/mosdepth

#! Make directories
## ${out_dir}/coverage/${build}/input must exist beforehand
mkdir -p ${outputDir}

#! Run commands

echo "Running mosdepth"

    mosdepth -t 4 -Q 20 -n -b "/path/to/retro/teDeletions.bed" \
    ${outputDir}/${sample_name}.${aligner}.del ${bamFile}

    mosdepth -t 4 -Q 20 -n -b "/path/to/retro/insertions.bed" \
    ${outputDir}/${sample_name}.${aligner}.ins ${bamFile}

echo "Finished mosdepth"

#! Deactivate ONT environment
conda deactivate
