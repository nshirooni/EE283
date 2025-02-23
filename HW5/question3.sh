srun -A CLASS-ECOEVO283 -c 1 --pty /bin/bash -i

module load vcftools/0.1.16 

vcftools --gzvcf output3.vcf.gz --chr X --from-bp 0 --to-bp 1000000 --012 --out snp_matrix

R

geno_data <- read.table("snp_matrix.012", header=FALSE)
geno_matrix <- as.matrix(geno_data[,-1])
geno_matrix <- t(geno_matrix)

png("snp_visualization.png", width=1200, height=800)
image(geno_matrix, 
      col=c("green", "yellow", "red"),
      breaks=c(-0.5, 0.5, 1.5, 2.5),
      xlab="SNP position",
      ylab="Samples",
      main="X Chromosome SNP Genotypes (First 1Mb)")

legend("topright", 
       legend=c("Ref/Ref", "Ref/Alt", "Alt/Alt"),
       fill=c("green", "yellow", "red"),
       border=NA)

dev.off()

