#! /bin/bash

survivor=/path/to/SURVIVOR

wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd152.GRCh38.variant_call.vcf.gz

gunzip nstd152.GRCh38.variant_call.vcf.gz

# We split the file per sample, get only insertions, and change chromosome numbers from "1" to "chr1"
truncate -s 0 nstd152.input
for i in \
    $(cat nstd152.GRCh38.variant_call.vcf | \
    awk '{match($8, "SAMPLE=([^;]+)", a); print a[1]}' | sort | uniq); do \
    
    cat nstd152.GRCh38.variant_call.vcf | grep -E "^#|SAMPLE=$i" | \
    grep -E "^#|SVTYPE=INS" | \
    awk -v'OFS=\t' \
    '{if ($0 ~ /^#/) {print $0} else { if ($1 ~ /^([[:digit:]]{1,2}|[XYM])/) {$1 = "chr"$1; print $0} } }' \
    > nstd152.$i.ins.vcf
    cat nstd152.$i.ins.vcf >> nstd152.input
done

$survivor merge nstd152.input 100 1 1 0 0 0 nstd152.ins.vcf

# We sortlike this because the header is not compatible with bcftools
cat \
    <(grep "^#" nstd152.ins.vcf) \
    <(cat nstd152.ins.vcf | grep -Ev "^#" | sort -k1V,1 -k2n,2) \
    > nstd152.ins.sorted.vcf
bgzip -f nstd152.ins.sorted.vcf > nstd152.ins.vcf.gz
tabix -p vcf --csi nstd152.ins.vcf.gz

python3 merge.py \
    -samples US HGSVC \
    -vcf all.me.vcf.gz nstd152.ins.vcf.gz -o us.vs.hgsvc.vcf \
    -d 60
