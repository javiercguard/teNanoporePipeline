#! /bin/bash

# requires bcftools

samples=$@

refGenome="GRCh38"
refInfix="hg38"

chromPattern=""

if [[ $refGenome = "GRCh38" ]] ; then
    chromPattern="^#|^chr([[:digit:]]+|[XY])	"
else
    chromPattern="^#|^([[:digit:]]+|[XY])	"
fi

survivor=/path/to/SURVIVOR
parentPath=/path/to/parentFolder
aligner="minimap"

for sample in ${samples[@]} ; do 

    sampleName=`echo $sample | cut -d"_" -f 1`

    echo $sampleName
    
    # SVIM + CuteSV
    # The better one is CuteSV so that one goes last

    truncate -s 0 survivorInput.tmp

    for file in "/path/to/${sample}/variant_calling/${refGenome}/"{svim,cutesv}"/minimap/"*${refInfix}*${aligner}*min3*gz ; do
        gunzip < $file > $(dirname $file)/$(basename $file ".gz")
        echo $(dirname $file)/$(basename $file ".gz") >> survivorInput.tmp
    done

    merged="$parentPath/$sample/variant_calling/$refGenome/$sampleName.$aligner.$refInfix.merged2.min3.vcf"
    $survivor merge survivorInput.tmp 10 1 1 0 0 1 ${merged}
    
    for file in "${parentPath}/${sample}/variant_calling/${refGenome}/"{svim,cutesv}"/minimap/"*${refInfix}*${aligner}*min3*gz ; do
        rm $(dirname $file)/$(basename $file ".gz")
    done
    rm survivorInput.tmp

    bcftools sort ${merged} -O v -o ${merged}

done
