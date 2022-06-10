source /path/to/conda.sh

conda activate pysam

for f in "folderName1" "folderName2" ; do

    sampleName=$(echo $1 | cut -f 1 -d_)

    python3 merge.py -samples cuteSV SVIM \
    -vcf /path/to/polish/$sampleName.hg38.minimap.{cutesv,svim}.min3.polish.vcf.gz \
    -d 20 \
    -o /path/to/data/$sampleName.merged.ins.min3.vcf

done

conda deactivate