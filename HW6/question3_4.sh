#Question 3 and 4
srun -A CLASS-ECOEVO283 -c 12 --pty /bin/bash -i 

mamba activate deeptools
module load samtools/1.15.1
module load bedtools2/2.30.0
module load ucsc-tools/v429

ref="ref/dm6.fa"
chromsizes="dm6.chrom.sizes"
Nreads=`samtools view -c -q 30 -F 4 condition1.sorted.bam`
Scale=`echo "1.0/($Nreads/1000000)" | bc -l`
samtools view -b condition1.sorted.bam | genomeCoverageBed -ibam - -g $ref -bg
-scale $Scale > condition1.coverage
various file types
bedGraphToBigWig condition1.coverage $chromsizes condition1.bw

samtools index condition1.sorted.bam 
bamCoverage -b condition1.sorted.bam -o condition1.bw --normalizeUsing RPKM --binSize 10 --extendReads 200 --ignoreDuplicates --maxFragmentLength 500

