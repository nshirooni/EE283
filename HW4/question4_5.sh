salloc -A ecoevo283 --ntasks=2 srun --pty /bin/bash -i

bamCoverage -b A4_filtered.bam -o A4_frag.bw --binSize 1 --normalizeUsing RPKM --extendReads
bamCoverage -b A5_filtered.bam -o A5_frag.bw --binSize 1 --normalizeUsing RPKM --extendReads

computeMatrix reference-point \
-S A4_frag.bw A5_frag.bw \
-R region.bed \
--referencePoint center \
-a 5000 -b 5000 \
-out matrix_frag.gz

bigwigCompare -b1 A4_frag.bw -b2 A5_frag.bw \
--operation ratio \
--region chrX:1904042:1904043 \
-o A4_vs_A5_frag_ratio.bw

