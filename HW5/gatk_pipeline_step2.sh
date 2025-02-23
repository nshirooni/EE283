#!/bin/bash
#SBATCH --job-name=hw5
#SBATCH -A CLASS-ECOEVO283
#SBATCH -p standard
#SBATCH --array=1-7  
#SBATCH --ntasks=1   
#SBATCH --cpus-per-task=4  
#SBATCH --mem=16G   

module load gatk/4.2.6.1
ref=/pub/nshiroon/EE283/DATA/ref/dmel-all-chromosome-r6.13.fasta
gvcf=/pub/nshiroon/EE283/DATA/DNAseq/allsample.g.vcf.gz  

mychr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" chrome.names.txt)

mkdir -p $HOME/scratch/$USER/tmp  

/opt/apps/gatk/4.2.6.1/gatk --java-options "-Xmx12g -Djava.io.tmpdir=$HOME/scratch/$USER/tmp" \
    GenotypeGVCFs \
    -R "$ref" \
    -V "$gvcf" \
    --intervals "$mychr" \
    -stand-call-conf 5 \
    -O "result.${mychr}.vcf.gz"

##Run outside of job
java -jar /opt/apps/picard-tools/2.27.1/picard.jar MergeVcfs \
$(printf 'INPUT=%s ' result.*.vcf.gz) \
O=all_variants.vcf.gz

