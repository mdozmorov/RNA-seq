#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=24
#PBS -M mdozmorov@vcu.edu
#PBS -m abe
#PBS -N STARindex
#PBS -j oe
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR

# Input folder
# DIRIN=01_trimmed_merged/
# Output folder
# DIROUT=02_subread-align

# Path to genome annotation files
## Human
DIRINDEX=/home/sequencing/data/ExtData/Ensembl/hg38
DNAINDEX=Homo_sapiens.GRCh38.dna.primary_assembly.fa
GENINDEX=Homo_sapiens.GRCh38.100.gtf
## Mouse
DIRINDEX=/home/mdozmorov/sequencing/data/ExtData/Ensembl/mm10
DNAINDEX=Mus_musculus.GRCm38.dna.primary_assembly.fa
GENINDEX=Mus_musculus.GRCm38.100.gtf

# mkdir -p $DIROUT

STAR --runThreadN 24 --runMode genomeGenerate --genomeDir $DIRINDEX --genomeFastaFiles $DIRINDEX"/"$DNAINDEX --sjdbGTFfile $DIRINDEX"/"$GENINDEX
