srun -A CLASS-ECOEVO283 -c 12 --pty /bin/bash -i 

module load subread/2.0.3
gtf="ref/dmel/dmel-all-r6.13.gtf"
myfile=`cat shortRNAseq.names.txt | tr "\n" " "`
featureCounts -p -T 12 -t exon -g gene_id -Q 30 -F GTF -a $gtf -o fly_counts.txt
$myfile 
