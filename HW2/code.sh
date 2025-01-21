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

mkdir /pub/nshiroon/EE283/HW2/RNAseq/
SourceDir="/data/class/ecoevo283/public/Bioinformatics_Course/RNAseq/RNAseq384plex_flowcell01"
DestDir="/pub/nshiroon/EE283/HW2/RNAseq/"
Labels="/data/class/ecoevo283/public/Bioinformatics_Course/RNAseq/RNAseq384_SampleCoding.txt"
tail -n +2 $Labels | head -n -1 > $DestDir/labels_rna.txt
File="/pub/nshiroon/EE283/HW2/RNAseq/labels_rna.txt"
while read -r number plex lane barcode plate row col well rest; do
    for plexdir in ${SourceDir}/Project_plex*/; do
        if [ -d "${plexdir}/Sample_${number}" ]; then
            READ1="${plexdir}/Sample_${number}/${number}_${plate}_${well}_R1.fastq.gz"
            READ2="${plexdir}/Sample_${number}/${number}_${plate}_${well}_R2.fastq.gz"
            ln -s "$READ1" "${DestDir}/$(basename $READ1)"
            ln -s "$READ2" "${DestDir}/$(basename $READ2)"
        fi
    done
done < "$File"
