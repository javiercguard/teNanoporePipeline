#!/bin/bash

#SBATCH -J SVIM
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
threads=16

if [ $build = "GRCh37" ]; then
 buildInfix="hg19"
 ref=$build_37
else
 buildInfix="hg38"
 ref=$build_38
fi

bamFile="${out_dir}/alignments/${build}/${aligner}/${sample_name}.${buildInfix}.bam"
vcfDir="${out_dir}/variant_calling/${build}/svim"
vcfTmp="${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.min1.unsorted.tmp.vcf"

#! Make directories
mkdir -p "${vcfDir}/${sample_name}_svim"
mkdir -p "${vcfDir}/${aligner}"

#! Run commands
echo "Starting SVIM at:"
date

conda activate "env2"

svim alignment --sequence_alleles --read_names "${vcfDir}"/"${sample_name}"_svim/ "${bamFile}" "${ref}"

conda deactivate

cp "${vcfDir}"/"${sample_name}"_svim/variants.vcf "${vcfTmp}"

conda activate ${env_path}

#! Sort the SVIM output
sortedVcf=${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.min1.sorted.tmp.vcf
bcftools sort ${vcfTmp} -O v -o ${sortedVcf}

#  Remove alleles for non-INS
cp ${sortedVcf} ${sortedVcf}.tmp
awk -v'OFS=\t' \
    '{if (substr($0, 1, 1) == "#") {print} else {match($8, /SVTYPE=([^;]+)/, a); if (a[1] == "INS") {print} else {$5="<"a[1]">";$4="N";print} } }' \
    ${sortedVcf}.tmp > ${sortedVcf}
rm ${sortedVcf}.tmp


# Filter by RE >= 2, and whatever support was listed in the config
# These come out sorted
for support in "3" "${sr_sniffles}" ; do 

	grep -E "^#|[	;]SUPPORT=([$(seq -s '' $support 9)]{1};?|[[:digit:]]{2,};?)" \
	${sortedVcf} \
	> "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.min${support}.sorted.tmp.vcf"

done

# Compress and index the vcfs, and remove the temporary ones
for support in "3" "${sr_sniffles}" ; do 

	definitiveVcf="${vcfDir}/${aligner}/${sample_name}.${buildInfix}.${aligner}.svim.min${support}.vcf.gz"

	bcftools view -O z \
	"${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.min${support}.sorted.tmp.vcf" \
	-o "${definitiveVcf}"
	bcftools index "${definitiveVcf}"

done

rm "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.svim.min"*".tmp.vcf"
rm -r "${vcfDir}"/"${sample_name}"_svim/

echo "Finishing SVIM at:"
date

conda deactivate