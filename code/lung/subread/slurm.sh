#!/bin/bash

#SBATCH --mem=32g 
#SBATCH --cpus-per-task=4
#SBATCH --time=7-00:00:00 
#SBATCH -p long
#SBATCH --begin=now+1hour

prj=/vast/projects/diffSplice/GibbsDTE-code

module load nextflow/23.04.2

nextflow pull plbaldoni/nextflow-rnaseq

nextflow run plbaldoni/nextflow-rnaseq -with-report report.html -with-trace -with-timeline timeline.html -resume -r main \
  --align --subjunc \
  --reads "$prj/data/lung/fastq/*{R1,R2}.fastq.gz" \
  --outdir "$prj/output/lung/subread/" \
  --subreadAnnoType "SAF" \
  --subreadAnno "$prj/data/annotation/hg38/230321-hg38_RefSeqStrict.saf.gz" \
  --subreadIndex "$prj/data/lung/index-subread/GRCh38.genome_index" \
  --subreadGenome "$prj/data/lung/index/GRCh38.p13.genome_sequins.fa.gz" \
  --gsize 3049315783

