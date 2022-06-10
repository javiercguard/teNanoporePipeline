source /path/to/conda.sh

conda activate pysam

python3 merge.py -samples US Indigen \
    -vcf /path/to/all.me.vcf.gz /path/to/nstd215.GRCh38.variant_call.prep.vcf.gz \
    -o /path/to/us.vs.indigen.vcf -d 60

python3 merge.py -samples US HGSVC \
    -vcf /path/to/all.me.vcf.gz /path/to/chaisson.ins.vcf.gz \
    -o /path/to/us.vs.hgsvc.vcf -d 60

python3 merge.py -samples US HGSVC Indigen \
    -vcf /path/to/all.me.vcf.gz /path/to/chaisson.ins.vcf.gz /path/to/nstd215.GRCh38.variant_call.prep.vcf.gz \
    -o /path/to/all3.vcf -d 60

conda deactivate