#!/bin/bash

#SBATCH --job-name=vcf_summary 
#SBATCH -A CLASS-ECOEVO283
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

source ~/.bashrc
mamba activate vcf_analysis

python3 vcf_summarize.py
