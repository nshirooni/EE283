module load fastqc/0.11.9
module loadtrimmomatic/0.39
java -jar /opt/apps/trimmomatic/0.39/trimmomatic-0.39.jar

cd /pub/nshiroon/EE283/HW2/DNAseq

SEQ1=/pub/nshiroon/EE283/HW2/DNAseq/ADL06_1_1.fq.gz
SEQ2=/pub/nshiroon/EE283/HW2/DNAseq/ADL06_1_2.fq.gz

fastqc $SEQ1 $SEQ2

java -jar /opt/apps/trimmomatic/0.39/trimmomatic-0.39.jar PE \
  -threads 12 \
  $SEQ1 $SEQ2 \
  trimmed_forward_paired.fastq.gz trimmed_forward_unpaired.fastq.gz \
  trimmed_reverse_paired.fastq.gz trimmed_reverse_unpaired.fastq.gz \
  ILLUMINACLIP:/opt/apps/trimmomatic/0.39/adapters/NexteraPE-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

SEQ1_trim=/pub/nshiroon/EE283/HW2/DNAseq/trimmed_forward_paired.fastq.gz
SEQ2_trim=/pub/nshiroon/EE283/HW2/DNAseq/trimmed_reverse_paired.fastq.gz

fastqc $SEQ1_trim $SEQ2_trim

cd /pub/nshiroon/EE283/HW2/ATACseq

ATAC_SEQ1=/pub/nshiroon/EE283/HW2/ATACseq/A4_ED_rep2_R1.fq.gz
ATAC_SEQ2=/pub/nshiroon/EE283/HW2/ATACseq/A4_ED_rep2_R1.fq.gz

fastqc $ATAC_SEQ1 $ATAC_SEQ2

java -jar /opt/apps/trimmomatic/0.39/trimmomatic-0.39.jar PE \
  -threads 12 \
  $ATAC_SEQ1 $ATAC_SEQ2 \
  trimmed_forward_paired.fastq.gz trimmed_forward_unpaired.fastq.gz \
  trimmed_reverse_paired.fastq.gz trimmed_reverse_unpaired.fastq.gz \
  ILLUMINACLIP:/opt/apps/trimmomatic/0.39/adapters/NexteraPE-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

ATAC_SEQ1_trim=/pub/nshiroon/EE283/HW2/ATACseq/trimmed_forward_paired.fastq.gz 
ATAC_SEQ2_trim=/pub/nshiroon/EE283/HW2/ATACseq/trimmed_reverse_paired.fastq.gz

mkdir /pub/nshiroon/EE283/HW2/ref
cd /pub/nshiroon/EE283/HW2/ref

ln -s /data/class/ecoevo283/public/Bioinformatics_Course/ref/dmel-all-chromosome-r6.13.fasta .
ln -s /data/class/ecoevo283/public/Bioinformatics_Course/ref/dmel-all-r6.13.gtf .

module load bwa/0.7.8
module load samtools/1.15.1
ref=/pub/nshiroon/EE283/HW2/ref/dmel-all-chromosome-r6.13.fasta

pathToRawData=/pub/nshiroon/EE283/HW3/data
pathToOutput=/pub/nshiroon/EE283/HW2/data/output

bwa index $ref

#Currently struggling with the 3rd step hoping to solve it soon!

