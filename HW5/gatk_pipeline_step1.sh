#!/bin/bash
#SBATCH --job-name=hw5
#SBATCH -A CLASS-ECOEVO283
#SBATCH -p standard
#SBATCH --array=1-4  
#SBATCH --ntasks=1    
#SBATCH --cpus-per-task=4  
#SBATCH --mem=16G    

module load java/1.8.0
module load gatk/4.2.6.1
module load picard-tools/2.27.1
module load samtools/1.15.1

ref=/pub/nshiroon/EE283/DATA/ref/dmel-all-chromosome-r6.13.fasta
cd /pub/nshiroon/EE283/DATA/DNAseq

samples=("ADL06" "ADL09" "ADL10" "ADL14")

sample=${samples[$SLURM_ARRAY_TASK_ID]}

java -jar /opt/apps/picard-tools/2.27.1/picard.jar AddOrReplaceReadGroups \
    I=${sample}.sort.bam \
    O=${sample}.RG.bam \
    SORT_ORDER=coordinate \
    RGPL=illumina \
    RGPU=D109LACXX \
    RGLB=Lib1 \
    RGID=${sample} \
    RGSM=${sample} \
    VALIDATION_STRINGENCY=LENIENT

java -jar /opt/apps/picard-tools/2.27.1/picard.jar MarkDuplicates \
    REMOVE_DUPLICATES=true \
    I=${sample}.RG.bam \
    O=${sample}.dedup.bam \
    M=${sample}_marked_dup_metrics.txt

samtools index ${sample}.dedup.bam

/opt/apps/gatk/4.2.6.1/gatk HaplotypeCaller \
    -R $ref \
    -I ${sample}.dedup.bam \
    --min-base-quality-score 30 \
    -ERC GVCF \
    -O ${sample}.g.vcf.gz

##Run outside of job
/opt/apps/gatk/4.2.6.1/gatk CombineGVCFs \
-R $ref $(for f in *.g.vcf.gz; do echo -V "$f"; done) \
-O allsample.g.vcf.gz

cat $ref | grep ">" | cut -f1 -d" " | tr -d ">" | head -n 7 >chrome.names.txt

