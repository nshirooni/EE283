#!/bin/bash
#SBATCH --job-name=hw2
#SBATCH -A CLASS-ECOEVO283
#SBATCH -p standard
#SBATCH --cpus-per-task=1


mkdir /pub/nshiroon/EE283/HW2/DNAseq/
SourceDir="/data/class/ecoevo283/public/Bioinformatics_Course/DNAseq/"
DestDir="/pub/nshiroon/EE283/HW2/DNAseq/"
FILES="$SourceDir/*"
for f in $FILES
do
    ff=$(basename $f)
    11
    echo "Processing $ff file..."
    ln -s $SourceDir/$ff $DestDir/$ff
done

mkdir /pub/nshiroon/EE283/HW2/ATACseq/
SourceDir="/data/class/ecoevo283/public/Bioinformatics_Course/ATACseq/"
DestDir="/pub/nshiroon/EE283/HW2/ATACseq/"
Labels="$SourceDir/README.ATACseq.txt"
head -n -3 $Labels | tail -n +2 > $DestDir/labels_for_atac.txt
File="/pub/nshiroon/EE283/HW2/ATACseq/labels_for_atac.txt"
while read p; do
    barcode=$(echo $p | cut -f1 -d" ")
    genotype=$(echo $p | cut -f2 -d" ")
    tissue=$(echo $p | cut -f3 -d" ")
    bioRep=$(echo $p | cut -f4 -d" ")
    
    READ1=$(find ${SourceDir}/ -type f -iname "*_${barcode}_R1.fq.gz")
    READ2=$(find ${SourceDir}/ -type f -iname "*_${barcode}_R2.fq.gz")
    
    ln -s "$READ1" "${DestDir}/${genotype}_${tissue}_rep${bioRep}_R1.fq.gz"
    ln -s "$READ2" "${DestDir}/${genotype}_${tissue}_rep${bioRep}_R2.fq.gz"
    
done < "$File"
