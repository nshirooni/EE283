#!/bin/bash
#SBATCH --job-name=hw1
#SBATCH -A CLASS-ECOEVO283
#SBATCH -p standard
#SBATCH --cpus-per-task=1

cd /pub/nshiroon/EE283/HW1
wget https://wfitch.bio.uci.edu/~tdlong/problem1.tar.gz
tar -xvf problem1.tar.gz
rm problem1.tar.gz

cd problem1

head -n 10 p.txt | tail -n 1
head -n 10 f.txt | tail -n 1

