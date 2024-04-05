#!/bin/bash
#SBATCH --mem=66g
#SBATCH --cpus-per-task=11
#SBATCH --time=12:00:00

module load salmon/1.10.2

salmonindex=../../../data/lung/index/transcripts_index/

for filename in ../../../data/lung/fastq/*_R1.fastq.gz; do

outname=${filename%_R1.fastq.gz}
outname=${outname##*/}

dt1=$(date '+%d/%m/%Y %H:%M:%S')

salmon quant \
-i $salmonindex \
-l A \
-1 $filename \
-2 ${filename/_R1.fastq.gz/_R2.fastq.gz} \
-p 10 \
--numBootstraps 100 \
--dumpEq \
-o ../../../output/lung/salmon-bootstrap/$outname # Real run
# -o /vast/scratch/users/baldoni.p/tmp/$outname # Timestamp run

dt2=$(date '+%d/%m/%Y %H:%M:%S')
echo "Timestamp Bootstrap Sample $outname: $dt1 - $dt2"

dt1=$(date '+%d/%m/%Y %H:%M:%S')

salmon quant \
-i $salmonindex \
-l A \
-1 $filename \
-2 ${filename/_R1.fastq.gz/_R2.fastq.gz} \
-p 10 \
--numGibbsSamples 100 \
--dumpEq \
-o ../../../output/lung/salmon-gibbs/$outname # Real run
# -o /vast/scratch/users/baldoni.p/tmp/$outname # Timestamp run

dt2=$(date '+%d/%m/%Y %H:%M:%S')
echo "Timestamp Gibbs Sample $outname: $dt1 - $dt2"

done
