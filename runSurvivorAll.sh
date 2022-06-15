#! /bin/bash

# requires bcftools

if [[ -z $1 ]] 
then 
    samples=("folderName1" "folderName2")
else 
    samples=$@
fi

refGenome="GRCh38"
refInfix="hg38"

chromPattern=""


# In hg38 chromosomes are called "chr1", but in hg19 it was "1"
if [[ $refGenome = "GRCh38" ]] ; then
    chromPattern="^#|^chr([[:digit:]]+|[XY])	"
else
    chromPattern="^#|^([[:digit:]]+|[XY])	"
fi

survivor=~/path/to/SURVIVOR
parentPath=/path/to/sampleFolders

aligner="minimap"

# echo "" > all.merged2.min3.input

truncate -s 0 survivorAll.tmp
truncate -s 0 survivorAllTrusty.tmp

for sample in ${samples[@]} ; do 

    sampleName=`echo $sample | cut -d"_" -f 1`

    echo $sampleName

    merged="$parentPath/$sample/variant_calling/$refGenome/$sampleName.$aligner.$refInfix.merged2.min3.vcf"

    echo $merged >> survivorAll.tmp

    mergedTrusty=$(dirname $merged)/$(basename $merged ".vcf").trusty.vcf
    cat $merged | grep -E "^#|SUPP=2;" > $mergedTrusty
    # cat $mergedTrusty | head -n600 | tail -n1
    echo $mergedTrusty >> survivorAllTrusty.tmp

done

echo "Merging all"

$survivor merge survivorAll.tmp 500 1 1 0 0 1 all.merged2.min3.vcf
bcftools sort all.merged2.min3.vcf -O v -o all.merged2.min3.vcf

echo "Merging trusty"

$survivor merge survivorAllTrusty.tmp 500 1 1 0 0 1 all.merged2.min3.trusty.vcf
bcftools sort all.merged2.min3.trusty.vcf -O v -o all.merged2.min3.trusty.vcf
