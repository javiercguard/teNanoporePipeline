#!/bin/bash

#SBATCH -J MINIMAP2
#SBATCH -A <<YOUR INFO>>
#SBATCH --nodes=3
#SBATCH --ntasks=48
#SBATCH --time=36:00:00
#SBATCH --mail-type=FAIL
#SBATCH --exclusive
#SBATCH -p <<YOUR INFO>>

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')

#! Modify the environment seen by the application. For this example we need the default modules.
. /etc/profile.d/modules.sh                # This line enables the module command
module purge                               # Removes all modules still loaded

#! Source config file
source $1

#! Set ONT enviroment
source $conda_path
ont_env=$env_path
conda activate $ont_env

#! Define variables
threads=16

if [ $build = "GRCh37" ]; then
 build_post="hg19"
 ref=$build_37
else
 build_post="hg38"
 ref=$build_38
fi

merged=${out_dir}/fastq/${sample_name}.merged.fastq
out_sam=${out_dir}/alignments/${build}/minimap/${sample_name}.${build_post}.sam
out_bam=${out_dir}/alignments/${build}/minimap/${sample_name}.${build_post}.bam

#! Make directories
mkdir -p ${out_dir}/alignments/${build}/minimap

#! Run commands
echo "Starting MINIMAP workflow at:"
date
echo "running minimap2"
minimap2 -x map-ont --MD -a -t $threads $ref $merged > $out_sam

samtools view -h $out_sam | samtools sort - -O BAM -o $out_bam

samtools index $out_bam

if [[ -f "$out_bam".bai ]]; then
  rm $out_sam
fi

echo "Finishing MINIMAP workflow at:"
date

#! Deactivate ONT environment
conda deactivate
