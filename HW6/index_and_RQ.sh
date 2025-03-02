#index and read quality
srun -A CLASS-ECOEVO283 -c 12 --pty /bin/bash -i 

module load samtools/1.15.1
ls *.chrX.bam | parallel -j 12 'samtools sort -@ 8 -o {.}.sorted.bam {} && samtools index {.}.sorted.bam'
ls *.sorted.chrX.bam | parallel -j 12 'samtools index {}'
ls *.sorted.chrX.bam | parallel -j 12 'samtools view -q 30 -b {} chrX | samtools sort -O BAM -o {.}.bam'
rm *.bai
ls *.sorted.chrX.bam | parallel -j 12 'samtools index {}'

