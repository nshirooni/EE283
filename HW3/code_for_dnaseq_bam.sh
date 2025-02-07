#!/bin/bash
#SBATCH --job-name=hw2
#SBATCH -A CLASS-ECOEVO283
#SBATCH -p standard
#SBATCH --cpus-per-task=2
#SBATCH --array=1-12

module load bwa/0.7.8
module load samtools/1.15.1

ref="/pub/nshiroon/EE283/HW2/ref/dmel-all-chromosome-r6.13.fasta"
file="dna_prefixes.txt"
prefix=`head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1`
pathToRawData="/pub/nshiroon/EE283/HW2/DNAseq"
pathToOutput="/pub/nshiroon/EE283/HW2/DNAseq"

bwa mem -t 2 -M $ref $pathToRawData/${prefix}_1.fq.gz $pathToRawData/${prefix}_2.fq.gz | 
samtools view -bS - > $pathToOutput/${prefix}.bam

samtools sort $pathToOutput/${prefix}.bam -o $pathToOutput/${prefix}.sort.bam

