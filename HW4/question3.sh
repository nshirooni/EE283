salloc -A ecoevo283 --ntasks=2 srun --pty /bin/bash -i


mamba create -n deeptools 

mamba activate deeptools

mamba install -c conda-forge -c bioconda deeptools

samtools index A4_filtered.bam
samtools index A5_filtered.bam

bamCoverage -b A4_filtered.bam -o A4.bw --binSize 1 --normalizeUsing RPKM
bamCoverage -b A5_filtered.bam -o A5.bw --binSize 1 --normalizeUsing RPKM

echo -e "X\t1904042\t1904043" > region.bed

computeMatrix reference-point \
-S A4.bw A5.bw \
-R region.bed \
--referencePoint center \
-a 5000 -b 5000 \
-out matrix.gz

plotProfile -m matrix.gz -out coverage_profile.png \
--plotTitle "Coverage near chrX:1,904,042"

