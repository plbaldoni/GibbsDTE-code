#!/bin/bash
#SBATCH --mem=66g
#SBATCH --cpus-per-task=11
#SBATCH --time=12:00:00

kallistoindex=../../../data/lung/index-kallisto/transcripts_index

mkdir ../../../output/lung/kallisto/

for filename in ../../../data/lung/fastq/*_R1.fastq.gz; do

outname=${filename%_R1.fastq.gz}
outname=${outname##*/}

dt1=$(date '+%d/%m/%Y %H:%M:%S')

~/lab_smyth/baldoni.p/software/kallisto/kallisto quant \
-i $kallistoindex \
--bootstrap-samples=100 \
--threads=10 \
-o ../../../output/lung/kallisto/$outname \
$filename ${filename/_R1.fastq.gz/_R2.fastq.gz}

dt2=$(date '+%d/%m/%Y %H:%M:%S')
echo "Timestamp Bootstrap Sample $outname: $dt1 - $dt2"

done
