#!/bin/bash

chr_ref="$1"
chr_query="$2"
query_gap_bed="$3"

module load SAMtools/1.9-foss-2016b
module load minimap2/2.17-foss-2018a

ref=$(basename $chr_ref |cut -f1 -d ".")
query=$(basename $chr_query |cut -f1 -d ".")

minimap2 -cx asm5 --cs -t24 $chr_ref $chr_query> ${chr_ref}_${chr_query}.paf  #self alignment by not setting primary chains and ignoring diagonal anchors -P -D
paftools.js liftover -q 0 -l 10000 ${chr_ref}_${chr_query}.paf $query_gap_bed > ${chr_ref}_${chr_query}.ngap.bed


