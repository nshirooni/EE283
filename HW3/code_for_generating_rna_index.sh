#!/bin/bash
#SBATCH --job-name=hw3
#SBATCH -A CLASS-ECOEVO283
#SBATCH -p standard
#SBATCH --cpus-per-task=8

module load hisat2/2.2.1
ref="/pub/nshiroon/EE283/DATA/ref/dmel-all-chromosome-r6.13.fasta"
gtf="/pub/nshiroon/EE283/DATA/ref/dmel-al-r6.13.gtf"
python hisat2_extract_splice_sites.py $gtf > dm6.ss
python hisat2_extract_exons.py $gtf > dm6.exon
hisat2-build -p 8 --exon dm6.exon --ss dm6.ss $ref dm6_trans
