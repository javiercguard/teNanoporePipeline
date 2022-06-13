#! /bin/bash

# Lets add a header and then compress + index

cat header.txt /path/to/all.me.txt > all.me.vcf

bcftools sort -O z -o /path/to/all.me.vcf.gz all.me.vcf
bcftools index /path/to/all.me.vcf.gz
