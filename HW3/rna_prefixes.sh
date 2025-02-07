#!/bin/bash
#SBATCH --job-name=hw3
#SBATCH -A CLASS-ECOEVO283
#SBATCH -p standard

samples="21148 21286 22162 21297 21029 22052 22031 21293 22378 22390"

for sample in $samples; do
    ls *${sample}*R1.fastq.gz | sed 's/_R1.fastq.gz//' >> rna_prefixes.txt
done

