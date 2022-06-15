#! /bin/bash

wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz

zcat hg38.fa.out.gz |  \
    tail -n+4 | \
    awk -v'OFS=\t' '{print $5,$6,$7,$9,$10,$11}' > repeatsReference.bed

zcat hg38.fa.out.gz | \
    tail -n+4 | \
    grep -E "(S|L)INE|Retroposon|LTR|DNA" | \
    awk -v'OFS=\t' '{print $5,$6,$7,$9,$10,$11}' >
    repeatsReferenceTE.bed

bgzip -k repeatsReferenceTE.bed
tabix -p bed --csi repeatsReferenceTE.bed.gz

cat repeatsReferenceTE.bed | \
    bedtools merge -c 4,5,6 -o collapse,collapse,collapse > 
    repeatsReferenceTENoOverlap.bed