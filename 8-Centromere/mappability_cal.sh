#!/bin/bash

bamfile="$1" #genomic Illumina mapping file, q20 as cutoff

cat $genome_fai |cut -f1,2 > chr.sizes.txt
bam_dir=$(dirname $bamfile)
cd $bam_dir
bam_prefix=$(basename $bamfile | cut -f1 -d ".")
# Reporting genome coverage for all positions in BEDGRAPH format.
module load BEDTools/2.28.0-foss-2018a
bedtools genomecov -ibam $bamfile -bga > ${bam_prefix}.Input.bedgraph
# Determine the read coverage cutoff by inspecting the bigwig format of bedgraph file [too low or too high]
module load ucsc/359
cat ${bam_prefix}.Input.bedgraph | sort -k1,1 -k2,2n > ${bam_prefix}.Input.sorted.bedgraph
# Extract regions of too low coverage (<x reads) and too high coverage (>y reads)
# 0-2/>100 reads covered region, with distance smaller than 1kb merged
cat ${bam_prefix}.Input.bedgraph |awk '{if($4<3||$4>100){print$0}}' |cut -f1,2,3 |sort -k1,1 -k2,2n |bedtools merge -i - -d 1000 > ${bam_prefix}_nonunique.bed
#cat ${bam_prefix}.Input.bedgraph |awk '{if($4<3||$4>100){print$0}}' |cut -f1,2,3 |sort -k1,1 -k2,2n |awk '{if($3-$2>=100){print$0}}' > ${bam_prefix}_nonunique.bed
