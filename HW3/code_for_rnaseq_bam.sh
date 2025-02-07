#!/bin/bash
#SBATCH --job-name=hw3
#SBATCH -A CLASS-ECOEVO283
#SBATCH -p standard
#SBATCH --cpus-per-task=8
#SBATCH --array=1-40
#SBATCH --mem=32G

module load samtools/1.15.1
module load hisat2/2.2.1

ref="/pub/nshiroon/EE283/DATA/ref/dm6_trans" 
file="/pub/nshiroon/EE283/HW3/rna_prefixes.txt" 
pathToRawData="/pub/nshiroon/EE283/DATA/RNAseq" 
pathToOutput="/pub/nshiroon/EE283/DATA/RNAseq/bam"
prefix=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $file)
R1="${pathToRawData}/${prefix}_R1.fastq.gz" 
R2="${pathToRawData}/${prefix}_R2.fastq.gz" 
unsorted_bam="${pathToOutput}/${prefix}.bam" 
sorted_bam="${pathToOutput}/${prefix}.sorted.bam"  
hisat2 -p 2 -x $ref -1 $R1 -2 $R2 | samtools view -bS - > $unsorted_bam
samtools sort $unsorted_bam -o $sorted_bam
samtools index $sorted_bam
rm $unsorted_bam
