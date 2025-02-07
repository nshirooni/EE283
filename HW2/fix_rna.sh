#!/bin/bash
#SBATCH --job-name=hw2
#SBATCH -A CLASS-ECOEVO283
#SBATCH -p standard
#SBATCH --cpus-per-task=1

SourceDir="/data/class/ecoevo283/public/Bioinformatics_Course/RNAseq/RNAseq384plex_flowcell01"
DestDir="/pub/nshiroon/EE283/DATA/RNAseq/"
Labels="/data/class/ecoevo283/public/Bioinformatics_Course/RNAseq/RNAseq384_SampleCoding.txt"
tail -n +2 $Labels | head -n -1 > $DestDir/labels_rna.txt
File="/pub/nshiroon/EE283/DATA/RNAseq/labels_rna.txt"


while read p;
do

        sampleno=$(echo $p | cut -f1 -d" ")
        i7index=$(echo $p | cut -f4 -d" ")
        lane=$(echo $p | cut -f3 -d" ")
        RILcode=$(echo $p | cut -f9 -d " ")
        TissueCode=$(echo $p | cut -f10 -d " ")
	Replicate=$(echo $p | cut -f11 -d " ")    
        
        READ1=$(find "${SourceDir}"/ -type f -iname "${sampleno}_${i7index}_${lane}_R1_001.fastq.gz")
        READ2=$(find "${SourceDir}"/ -type f -iname "${sampleno}_${i7index}_${lane}_R2_001.fastq.gz")
         
        ln -s "$READ1" "${DestDir}/${RILcode}_${TissueCode}_${Replicate}_R1.fastq.gz"
        ln -s "$READ2" "${DestDir}/${RILcode}_${TissueCode}_${Replicate}_R2.fastq.gz"

        
done < $File

