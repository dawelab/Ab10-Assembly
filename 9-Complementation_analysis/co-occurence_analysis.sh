#!/bin/bash

#bed files need to be sorted by sort -k1,1 -k2,2n
genome="$1"
bed1="$2"
bed2="$3"

module load SAMtools/1.9-foss-2016b
samtools faidx $genome 
genome_prefix=$(basename $genome |cut -f1 -d ".")
cat ${genome_prefix}.fai |cut -f1,2 > ${genome_prefix}.genome

module load BEDTools/2.28.0-foss-2018a
'''Co-occurence analyses'''
cat $bed1 | bedtools fisher -a - -b $bed2 -g ${genome_prefix}.genome

