#!/bin/bash

#SBATCH -J NGMLR
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

fastq=${out_dir}/fastq/${sample_name}.fastq.fof
merged=${out_dir}/fastq/${sample_name}.merged.fastq
out_sam=${out_dir}/alignments/${build}/ngmlr/${sample_name}.${build_post}.sam
out_bam=${out_dir}/alignments/${build}/ngmlr/${sample_name}.${build_post}.bam

#! Make directories 
mkdir -p ${out_dir}/alignments/${build}/ngmlr/

#! Run commands
echo "Starting NGMLR workflow at:"
date

# echo "mapping reads with ngmlr"
ngmlr -t $threads -r $ref -q $merged -o $out_sam -x ont

echo "converting sam to bam"
samtools view -h $out_sam | samtools sort - -O BAM -o $out_bam

samtools index $out_bam

if [[ -f "$out_bam" ]]; then
  rm $out_sam
fi

echo "Finishing NGMLR workflow at:"
date

#! Deactivate ONT encironment
conda deactivate
