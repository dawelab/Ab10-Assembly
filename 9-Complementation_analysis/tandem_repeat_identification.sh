#!/bin/bash

genome="$1"

module load SAMtools/1.9-foss-2016b
module load minimap2/2.17-foss-2018a

# chr level self alignment 
samtools faidx $genome
for chr in chr{1..10}; do 
samtools faidx $genome $chr > ${chr}.fa #extract chr fa file 
minimap2 -t24 -PD -k19 -w19 -m200 ${chr}.fa ${chr}.fa > ${chr}_self.paf #self alignment by not setting primary chains and ignoring diagonal anchors -P -D
done

module load BEDTools/2.26.0-foss-2016b

# Extract tandem repeat from self alignment paf files (>25 Kb repetitiveness)
for file in *_self.paf;do
name=$(basename $file |cut -f1 -d ".")
chr=$(basename $file |cut -f1 -d "." | cut -f1 -d "_")
cat $file | awk '{if($11>25000){print$0}}' > ${name}.filtered.paf
cat ${name}.filtered.paf|sort -k3,3n |cut -f1,3,4 |bedtools merge -i - -d 25000 > ${chr}.bed
done
cat *.bed > v2.repeat2.bed

# Repeat array 
cat v2.repeat2.bed |sort -k1,1 -k2,2n|bedtools merge -i - -d 300000 > repeat.array.bed

