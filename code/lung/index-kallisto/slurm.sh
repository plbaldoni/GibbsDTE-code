#!/bin/bash

#SBATCH --mem=64g
#SBATCH --time=24:00:00

# Creating kallisto index for hg38

cp ../../../data/lung/index/gencode.v33.transcripts_sequins.fa.gz ./

gunzip gencode.v33.transcripts_sequins.fa.gz

~/lab_smyth/baldoni.p/software/kallisto/kallisto index \
-i ../../../data/lung/index-kallisto/transcripts_index \
gencode.v33.transcripts_sequins.fa

# Removing extra files

rm -rf gencode.v33.transcripts_sequins.fa
