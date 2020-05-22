#!/bin/bash

module load BEDTools/2.28.0-foss-2018a
#all three bed files need to be sorted by sort -k1,1 -k2,2n
bed1="$1"
bed2="$2" 
bed3="$3"

# in this case, bed1 is ngap bed file (13n, 100n, >10kb), bed2 is repeat array / heterozygosity bed file
'''map Ngap to tandem repeat array and heterozygosity'''
#Extract fully assembled repeats, ont-compelmented repeats, and ont-failed repeats'''
windowBed -a $bed1 -b $bed2 -r 500000 -l 500000 | cut -f4,5,6 |sort -k1,1 -k2,2n |uniq > unassembled.bed
bedtools intersect -a $bed2 -b unassembled.bed -v |sort -k1,1 -k2,2n |uniq > assembled.bed
bedtools intersect -a unassembled.bed -b $bed3 -u |sort -k1,1 -k2,2n |uniq > uncomplemented.bed
bedtools intersect -a unassembled.bed -b uncomplemented.bed -v |sort -k1,1 -k2,2n |uniq > complemented.bed

'''Quantify the distribution of 100n, 13n and >=10000 overlapped with tandem repeats'''
windowBed -a $bed1 -b $bed2 -l 500000 -r 500000 |cut -f1,2,3 |sort -k2,2n |uniq |awk '{if($3-$2>10000){print$0}}' |wc -l
windowBed -a $bed1 -b $bed2 -l 500000 -r 500000 |cut -f1,2,3 |sort -k2,2n |uniq |awk '{if($3-$2=="13"){print$0}}' |wc -l
windowBed -a $bed1 -b $bed2 -l 500000 -r 500000 |cut -f1,2,3 |sort -k2,2n |uniq |awk '{if($3-$2=="100"){print$0}}' |wc -l

