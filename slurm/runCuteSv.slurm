#!/bin/bash

#SBATCH -J CUTESV
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

#! Enable conda usage, particularly my CuteSV env
source $conda_path
conda activate "env2"

#! Define variables
aligner=$2
threads=16

if [ $build = "GRCh37" ]; then
 buildInfix="hg19"
 ref=$build_37
else
 buildInfix="hg38"
 ref=$build_38
fi

bamFile="${out_dir}/alignments/${build}/${aligner}/${sample_name}.${buildInfix}.bam"
vcfDir="${out_dir}/variant_calling/${build}/cutesv/${aligner}"
vcfTmpSeq="${vcfDir}/${sample_name}.${buildInfix}.${aligner}.cutesv.min1.unsorted.tmp.withSeq.vcf"
vcfTmp="${vcfDir}/${sample_name}.${buildInfix}.${aligner}.cutesv.min1.unsorted.tmp.vcf"

#! Make directories
mkdir -p "${vcfDir}"

#! Run commands
echo "Starting CuteSV at:"
date

cuteSV -t ${threads} --max_cluster_bias_INS 100 -s2 -L -1 \
    --report_readid \
    --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
    ${bamFile} ${ref} ${vcfTmpSeq} ${vcfDir}

conda deactivate

conda activate ${env_path}

awk -v'OFS=\t' \
    '{if (substr($0, 1, 1) == "#") {print} else {match($8, /SVTYPE=([^;]+)/, a); if (a[1] == "INS") {print} else {$5="<"a[1]">";$4="N";print} } }' \
    ${vcfTmpSeq} > ${vcfTmp}
rm ${vcfTmpSeq}

#! Sort the output
sortedVcf=${vcfDir}/${sample_name}.${buildInfix}.${aligner}.cutesv.min1.sorted.tmp.vcf
bcftools sort "$vcfTmp" -O v -o "${sortedVcf}"

#! Filter by RE

for support in "2" "${sr_sniffles}" ; do 
    grep -E "^#|[	;]RE=([$(seq -s '' $support 9)]{1};?|[[:digit:]]{2,};?)" \
    ${sortedVcf} \
    > "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.cutesv.min${support}.sorted.tmp.vcf"
done

# Compress and index the vcfs, and remove the temporary ones
for support in "1" "2" "${sr_sniffles}" ; do 

    definitiveVcf="${vcfDir}/${sample_name}.${buildInfix}.${aligner}.cutesv.min${support}.vcf.gz"

    bcftools view -O z \
    "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.cutesv.min${support}.sorted.tmp.vcf" \
    -o "${definitiveVcf}"
    bcftools index "${definitiveVcf}"

    if [[ -f "${definitiveVcf}" ]]; then
        rm "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.cutesv.min${support}."*".tmp.vcf"
    fi

done
echo "Finishing CuteSV at:"
date

#! Deactivate ONT environment
conda deactivate