source /path/to/conda.sh

conda activate env3

python3 ~/datos-javier/ins/merge.py -samples \
    sample1 sample2 \
    -vcf /path/to/sample1.merged.ins.min3.vcf.gz /path/to/sample1.merged.ins.min3.vcf.gz \
    -o /path/to/all.merged.ins.85.min3.vcf \
    -d 60 --consider-sequence

conda deactivate

conda activate env4 

python3 checkInsertionsMultiSample.py /path/to/all.merged.ins.85.min3.vcf .

conda deactivate
