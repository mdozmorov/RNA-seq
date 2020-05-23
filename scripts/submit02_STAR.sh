#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=4
#PBS -M mdozmorov@vcu.edu
#PBS -m abe
#PBS -N STAR
#PBS -j oe
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR

# Input folder
DIRIN=01_trimmed/
# Output folder
DIROUT=02_STAR-align
# Path to genome annotation files
## Human
# DIRINDEX=/home/sequencing/data/ExtData/Ensembl/hg38
# DNAINDEX=Homo_sapiens.GRCh38.dna.primary_assembly.fa
# GENINDEX=Homo_sapiens.GRCh38.100.gtf
## Mouse
DIRINDEX=/home/mdozmorov/sequencing/data/ExtData/Ensembl/mm10
DNAINDEX=Mus_musculus.GRCm38.dna.primary_assembly.fa
GENINDEX=Mus_musculus.GRCm38.100.gtf

mkdir -p $DIROUT

# List of input files, each string contains (comma-separated) file name(s)
# Space separates first and second read pairs
# find 01_trimmed/ -type f -name "*val_1.fq.gz" | sort > input01_toStarAlign.list
# find 01_trimmed/ -type f -name "*val_2.fq.gz" | sort >> input01_toStarAlign.list
INPUT=input01_toStarAlign.list


# Paired end
echo "Starting STAR alignment on "`date`

cat $INPUT | while read read1 read2
do
	prefix=$(echo $read1 | cut -f1 -d "," | xargs basename | sed 's/_L006_R1_001_val_1\.fq\.gz//')
	echo "Starting alignment of $prefix on "`date` 
	STAR --runThreadN 4 --genomeDir $DIRINDEX --readFilesIn $read1 $read2 --readFilesCommand zcat --outSAMtype BAM Unsorted --outSAMorder Paired --outReadsUnmapped Fastx --outFileNamePrefix $DIROUT/$prefix. --sjdbGTFfile $DIRINDEX"/"$GENINDEX --quantMode GeneCounts --outFilterMultimapNmax 10 --bamRemoveDuplicatesType UniqueIdenticalNotMulti
	echo "Finished alignment of $prefix on "`date`

done

echo "Completed STAR alignment on "`date`

# Duplicate removal: https://github.com/alexdobin/STAR/issues/175
