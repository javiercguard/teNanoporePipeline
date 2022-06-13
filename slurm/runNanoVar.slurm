#!/bin/bash

#SBATCH -J NANOVAR
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

#! Enable conda usage
source $conda_path

#! Define variables
aligner=$2

if [ $build = "GRCh37" ]; then
 buildInfix="hg19"
 ref=$build_37
else
 buildInfix="hg38"
 ref=$build_38
fi

bamFile="${out_dir}/alignments/${build}/${aligner}/${sample_name}.${buildInfix}.q20.bam"
vcfDir="${out_dir}/variant_calling/${build}/nanovar/${aligner}"
tmpDir="${vcfDir}/tmp"
vcfTmp="${vcfDir}/${sample_name}.${buildInfix}.${aligner}.nanovar.unsorted.tmp.vcf"

conda activate ${env_path}

fullBamFile="${out_dir}/alignments/${build}/${aligner}/${sample_name}.${buildInfix}.bam"
samtools view -hb -q 20 ${fullBamFile} > ${bamFile}
samtools index ${bamFile}

conda deactivate

#! Make directories
mkdir -p "${tmpDir}"

#! Run commands

conda activate "/rds/user/bd417/hpc-work/conda_envs/nanovar"

echo "Starting NanoVar at:"
date

nanovar -t 15 -c 1 "${bamFile}" "${ref}" "${tmpDir}"

conda deactivate

rm "${bamFile}"
rm "${bamFile}.bai"

conda activate ${env_path}

mv "${tmpDir}/${sample_name}.${buildInfix}.q20.nanovar.pass.vcf" ${vcfTmp}

#! Sort the NanoVar output

sortedVcf="${vcfDir}/${sample_name}.${buildInfix}.${aligner}.nanovar.min1.sorted.tmp.vcf"
bcftools sort "$vcfTmp" -O v \
-o "${sortedVcf}"

#! Filter by SR

for support in "2" "${sr_sniffles}" ; do 
    grep -E "^#|[	;]SR=([$(seq -s '' $support 9)]{1};?|[[:digit:]]{2,};?)" \
    ${sortedVcf} \
    > "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.nanovar.min${support}.sorted.tmp.vcf"
done

# Compress and index the vcfs, and remove the temporary ones
for support in "1" "2" "${sr_sniffles}" ; do 

    definiteDir="${out_dir}/variant_calling/${build}/nanovar/${aligner}"
    definitiveVcf="${definiteDir}/${sample_name}.${buildInfix}.${aligner}.nanovar.min${support}.vcf.gz"

    bcftools view -O z \
    "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.nanovar.min${support}.sorted.tmp.vcf" \
    -o "${definitiveVcf}"
    bcftools index "${definitiveVcf}"

done

#! Compress and index the vcfs, and remove the temporary ones. The naming convention make looping pointless

echo "Finishing NanoVar at:"
date

#! Deactivate ONT environment
conda deactivate