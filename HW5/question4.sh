srun -A CLASS-ECOEVO283 -c 1 --pty /bin/bash -i

awk '{
     zero=0;
     one=0;
     two=0;
     for(i=2; i<=NF; i++) {
         if($i == 0) zero++;
         else if($i == 1) one++;
         else if($i == 2) two++;
     }
     if(zero >= 2 && two >= 2) 
     
R

geno_data <- read.table("filtered_snp_matrix.012", header=FALSE)
geno_matrix <- as.matrix(geno_data[,-1])

geno_matrix <- t(geno_matrix)

png("filtered_snp_visualization.png", width=1200, height=800)
image(geno_matrix, 
      col=c("green", "red"),  # Only need green and red now
      breaks=c(-0.5, 1.5, 2.5),  # Adjusted breaks for just 0 and 2
      xlab="SNP position",
      ylab="Samples",
      main="X Chromosome SNP Genotypes - Filtered for 0/0 vs 1/1 pattern")

legend("topright", 
       legend=c("Ref/Ref (0/0)", "Alt/Alt (1/1)"),
       fill=c("green", "red"),
       border=NA)

dev.off()
