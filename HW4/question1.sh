#!/bin/bash
#SBATCH --job-name=hw3
#SBATCH -A CLASS-ECOEVO283
#SBATCH -p standard

module load samtools/1.15.1
module load bedtools2/2.30.0

cd /pub/nshiroon/EE283/DATA/DNAseq

A4="/pub/nshiroon/EE283/DATA/DNAseq/ADL06.bam"
A5="/pub/nshiroon/EE283/DATA/DNAseq/ADL09.bam"

samtools sort -o ADL06_sorted.bam $A4
samtools sort -o ADL09_sorted.bam $A5

A4="/pub/nshiroon/EE283/DATA/DNAseq/ADL06_sorted.bam"
A5="/pub/nshiroon/EE283/DATA/DNAseq/ADL09_sorted.bam"

samtools index $A4
samtools index $A5

samtools view -b -q 30 $A4 X:1880000-2000000 > A4_filtered.bam
samtools view -b -q 30 $A5 X:1880000-2000000 > A5_filtered.bam
