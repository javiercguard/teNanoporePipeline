#!/bin/bash

#SBATCH -J LRA
#SBATCH -A <<YOUR ACCOUNT>>
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=36:00:00
#SBATCH --mail-type=FAIL
#SBATCH --exclusive
#SBATCH -p <<YOUR CLUSTER>>

#! Set up your modules, if any
. /etc/profile.d/modules.sh                # This line enables the module command
module purge                               # Removes all modules still loaded

#! Source config file
source $1

#! Set ONT enviroment 
source $conda_path
conda activate "envLra"

#! Define variables
threads=63

if [ $build = "GRCh37" ]; then
 build_post="hg19"
 ref=$build_37
else
 build_post="hg38"
 ref=$build_38
fi

fastq=${out_dir}/fastq/${sample_name}.fastq.fof
merged=${out_dir}/fastq/${sample_name}.merged.fastq
out_sam=${out_dir}/alignments/${build}/lra/${sample_name}.${build_post}.sam
out_bam=${out_dir}/alignments/${build}/lra/${sample_name}.${build_post}.bam

#! Make directories 
mkdir -p ${out_dir}/alignments/${build}/lra/

#! Run commands
echo "Starting LRA workflow at:"
date

lra align -ONT --printMD -t ${threads} -p s ${ref} ${merged} > ${out_sam}

echo "converting sam to bam"
samtools view -h $out_sam | samtools sort - -O BAM -o $out_bam

samtools index $out_bam

if [[ -f "$out_bam" ]]; then
  rm $out_sam
fi

echo "Finishing LRA workflow at:"
date

#! Deactivate ONT encironment
conda deactivate
