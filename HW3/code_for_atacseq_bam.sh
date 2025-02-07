#!/bin/bash
#SBATCH --job-name=hw3
#SBATCH -A CLASS-ECOEVO283
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --array=1-24

module load bwa/0.7.8
module load samtools/1.15.1

ref="/pub/nshiroon/EE283/DATA/ref/dmel-all-chromosome-r6.13.fasta"
file="atac_prefixes.txt"
prefix=`head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1`
pathToRawData="/pub/nshiroon/EE283/DATA/ATACseq"
pathToOutput="/pub/nshiroon/EE283/DATA/ATACseq"

bwa mem -t 2 -M $ref $pathToRawData/${prefix}_R1.fq.gz $pathToRawData/${prefix}_R2.fq.gz |
samtools view -bS - > $pathToOutput/${prefix}.bam

samtools sort $pathToOutput/${prefix}.bam -o $pathToOutput/${prefix}.sort.bam
