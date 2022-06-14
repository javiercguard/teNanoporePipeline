# 
# This is an example configuration file similar to the ones taht were used for the publication.
# 

#[DATA SAMPLE]
#Full path to the fastq files directoriess. If multiple, comma sepparated
data_dirs=/path/to/fastq_pass/,/path/to/fastq_pass2/,/path/to/fastq_fail/

#Full path to the output directory
out_dir=/path/to/outpurDirectory/
#Sample name
sample_name=exampleName

#[ANALYSIS]
#Genome build. Options: GRCh37,GRCh38
build=GRCh38
#Read evidence option, not only for Sniffles
sr_sniffles=5 

#[MISCELLANEOUS]
#Full path to condaSetUp.sh
conda_path=/path/to/condaSetUp.sh
#Full path to environment
env_path=/path/to/envs/ont
#Full path to the scripts directory
scripts_dir=/path/to/slurm_commands # ../slurmScripts in the presented setup
#Full path to the b37 reference genome
build_37=/path/to/Homo_sapiens.GRCh37.75.dna.fasta
#Full path to the b38 reference genome
build_38=/path/to/Homo_sapiens.GRCh38.dna.fasta
#Full path to the b37 bedfile
ref_bed_37=/path/to/human_hg19.bed
#Full path to the b38 bedfile
ref_bed_38=/path/to/human_hg38.bed