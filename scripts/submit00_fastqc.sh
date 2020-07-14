#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -l nodes=1:ppn=12
#PBS -M mdozmorov@vcu.edu
#PBS -m abe
#PBS -N fastqc
#PBS -j oe
# PBS -o /path/to/stderr-stdout/output

cd $PBS_O_WORKDIR

# Input folder
DIRIN=00_raw/

# Output folder
DIROUT=00_fastqc-raw

mkdir -p $DIROUT

for file in `find $DIRIN -type f -name "*.gz" | sort`; do 
	fastqc -t 12 -o $DIROUT --noextract $file; 
done

# Extract tab-separated FASTQC summary
# python fastqc-summary -s $DIROUT > $DIROUT".txt"