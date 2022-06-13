#!/bin/bash
#!
## author: Belen de la Morena Barrio bd417@cam.ac.uk
## updated last: 16 Jul 2019

#SBATCH -J Sniffles
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
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

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
vcfDir="${out_dir}/variant_calling/${build}/sniffles"

echo ${vcfDir}

vcfTmp="${vcfDir}/${sample_name}.${buildInfix}.${aligner}.sniffles.min2.unsorted.tmp.vcf"

#! Make directories
mkdir -p "${out_dir}/variant_calling/${build}/sniffles"

#! Run commands
echo "Starting SNIFFLES at:"
date


/path/to/Sniffles-master/bin/sniffles-core-1.0.12/sniffles -t ${threads} -s 2 -m ${bamFile} -v ${vcfTmp}

# Sort the Sniffles output

conda activate "env1"

sortedVcf=${vcfDir}/${sample_name}.${buildInfix}.${aligner}.sniffles.min2.sorted.tmp.vcf
bcftools sort ${vcfTmp} -O v -o ${sortedVcf}

# Filter by RE >= 2, and whatever support was listed in the config
# These come out sorted
for support in "3" "${sr_sniffles}" ; do 

	grep -E "^#|SUPPORT=([$(seq -s '' $support 9)]{1};?|[[:digit:]]{2,};?)" \
	${sortedVcf} \
	> "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.sniffles.min${support}.sorted.tmp.vcf"

done

# Compress and index the vcfs, and remove the temporary ones
for support in "2" "3" "${sr_sniffles}" ; do 

	definitiveVcf="${vcfDir}/${sample_name}.${buildInfix}.${aligner}.sniffles.min${support}.vcf.gz"

	bcftools view -O z \
	"${vcfDir}/${sample_name}.${buildInfix}.${aligner}.sniffles.min${support}.sorted.tmp.vcf" \
	-o "${definitiveVcf}"
	bcftools index "${definitiveVcf}"

	if [[ -f "${definitiveVcf}" ]]; then
		rm "${vcfDir}/${sample_name}.${buildInfix}.${aligner}.sniffles.min${support}."*".tmp.vcf"
	fi

done

echo "Finishing SNIFFLES at:"
date

#! Deactivate ONT environment
conda deactivate
