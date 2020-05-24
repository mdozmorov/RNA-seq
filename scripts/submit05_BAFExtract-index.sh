#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -M mdozmorov$vcu.edu
#PBS -m abe
#PBS -N BAFindex
#PBS -j oe
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR


# Path to ENSEMBL genome annotation files
## Human
# DIRINDEX=/home/sequencing/data/ExtData/Ensembl/hg38
# DNAINDEX=Homo_sapiens.GRCh38.dna.primary_assembly.fa
# GENINDEX=Homo_sapiens.GRCh38.100.gtf
# ## Mouse
# # DIRINDEX=/home/mdozmorov/sequencing/data/ExtData/Ensembl/mm10
# # DNAINDEX=Mus_musculus.GRCm38.dna.primary_assembly.fa
# # GENINDEX=Mus_musculus.GRCm38.100.gtf

# genome_fasta_pileup_dir=/home/sequencing/data/ExtData/BAFExtract/hg38/
# BAFExtract -preprocess_FASTA $DIRINDEX"/"$DNAINDEX $genome_fasta_pileup_dir

# Path to UCSC genome annotation files
# DIRINDEX=/home/sequencing/data/ExtData/UCSC/mm10
DIRINDEX=/home/mdozmorov/sequencing/data/ExtData/UCSC/hg38gdc
DNAINDEX=GRCh38.d1.vd1.fa

genome_fasta_pileup_dir=/home/sequencing/data/ExtData/BAFExtract/hg38gdc/
BAFExtract -preprocess_FASTA $DIRINDEX"/"$DNAINDEX $genome_fasta_pileup_dir
