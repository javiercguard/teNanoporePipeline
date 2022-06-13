#! /bin/bash

wget https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/nstd215.GRCh38.variant_call.vcf.gz

# We need to change the CHROM column from "1" to "chr1"
cat nstd215.GRCh38.variant_call.vcf | awk -v'OFS=\t' \
    '{if ($0 ~ /^#/) {print $0} else { if ($1 ~ /^([[:digit:]]{1,2}|[XYM])/) {$1 = "chr"$1; print $0} } }' | \
    awk  -v'OFS=\t' \
    '{if ($0 ~ /^#/) {print $0} else {$5="<INS>"; print $0}}' > \
    nstd215.GRCh38.variant_call.prep.vcf

bgzip -fk nstd215.GRCh38.variant_call.prep.vcf
tabix -f --csi nstd215.GRCh38.variant_call.prep.vcf.gz

# python3 merge.py \
#     -samples US Indigen \
#     -vcf nstd215.GRCh38.variant_call.prep.vcf.gz -o us.vs.hgsvc.vcf \
#     -d 60
