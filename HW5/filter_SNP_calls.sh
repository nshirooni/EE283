srun -A CLASS-ECOEVO283 -c 1 --pty /bin/bash -i 

module load bcftools/1.15.1

bcftools filter -i 'FS<40.0 && SOR<3 && MQ>40.0 && MQRankSum>-5.0 && MQRankSum<5 && QD>2.0 && ReadPosRankSum>-4.0 && INFO/DP<16000' \
-O z -o output1.vcf.gz all_variants.vcf.gz 

bcftools filter -S . -e 'FMT/DP<3 | FMT/GQ<20' \
-O z -o output2.vcf.gz output1.vcf.gz 

bcftools filter -e 'AC==0 || AC==AN' --SnpGap 10 output2.vcf.gz | \
bcftools view -m2 -M2 -v snps -O z -o output3.vcf.gz

