#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -M mdozmorov$vcu.edu
#PBS -m abe
#PBS -N BAFExtract
#PBS -j oe
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR

# Input folder
DIRIN=02_subread-align/

# Output folder
DIROUT=05_BAFExtract

mkdir -p $DIROUT

genome_list=/home/sequencing/data/ExtData/BAFExtract/hg38.list
genome_fasta_pileup_dir=/home/sequencing/data/ExtData/BAFExtract/hg38/

for file in `find $DIRIN -type f -name "*.bam"`; do
	# echo $file;
	bam_file=$file;
	output_baf_file=$DIROUT"/"`basename $file .bam`".snp";
	samtools view $bam_file | BAFExtract -generate_compressed_pileup_per_SAM stdin $genome_list $DIROUT 50 0;
	BAFExtract -get_SNVs_per_pileup  $genome_list $DIROUT $genome_fasta_pileup_dir 20 4 0.1 $output_baf_file;
done


