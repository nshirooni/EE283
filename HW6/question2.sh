#Question 2
srun -A CLASS-ECOEVO283 -c 12 --pty /bin/bash -i 

module load samtools/1.15.1
samtools merge -o condition1.bam A4_WD_rep1.sorted.bam A4_WD_rep2.sorted.bam A4_WD_rep4.sorted.bam
module load bedtools2/2.30.0
bedtools bamtobed -i condition1.bam | \
  awk -F $'\t' 'BEGIN {OFS = FS} {
      if ($6 == "+") {
          $2 = $2 + 4  # Adjust start for + strand
      } else if ($6 == "-") {
          $3 = $3 - 5  # Adjust end for - strand
      }
      print $0
  }' > condition1.tn5.bed 

module load macs/2.2.7.1
macs2 callpeak -t condition1.tn5.bed -n condition1 -f BED -g mm -q 0.01 \
  --nomodel --shift -75 --extsize 150 --call-summits --keep-dup all -B

macs2 callpeak -t condition1.tn5.bed -n condition1 -f BED -g mm -q 0.01 \
  --nomodel --shift -75 --extsize 150 --keep-dup all -B --broad

module load ucsc-tools/v429

rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64.v369/fetchChromSizes .
./fetchChromSizes dm6 > dm6.chrom.sizes

LC_COLLATE=C sort -k1,1 -k2,2n ?? .broad_treat_pileup.bdg > ?? .broad_treat_pileup.sorted.bdg

bedGraphToBigWig ?? .broad_treat_pileup.sorted.bdg dm6.chrom.sizes ?? .broad_peaks.b

