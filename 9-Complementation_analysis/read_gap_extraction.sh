#!/bin/bash

module load BEDTools/2.28.0-foss-2018a
genome="$1"
reads="$2" 

reads_prefix=$(echo $reads |cut -f1 -d ".")
genome_prefix=$(basename $genome |cut -f1 -d ".")


module load seqtk
seqtk seq -L 10000 $reads |gzip > ${reads_prefix}.10kb.fa.gz

ml minimap2/2.17-foss-2018a
#Align all reads to ref genome
minimap2 -ax map-pb $genome ${reads_prefix}.10kb.fa.gz > ${genome_prefix}.sam 
samtools view -b -@ 12 -o ${genome_prefix}.bam ${genome_prefix}.sam 
samtools sort -@ 12 -o ${genome_prefix}.sorted.bam ${genome_prefix}.bam 
samtools index ${genome_prefix}.sorted.bam

#Read mapping inspection / display
module load Python/3.6.4-foss-2018a
module load deepTools/3.2.1-foss-2018a-Python-3.6.4
bamCoverage --bam ${genome_prefix}.sorted.bam --binSize 10000 --numberOfProcessors max --normalizeUsing CPM --ignoreDuplicates --minMappingQuality 0 --outFileFormat bedgraph --outFileName ${genome_prefix}.10kb.bedgraph

'''obtain read gap'''
# Identify read distribution by checking unaligned positions 
module load BEDTools/2.28.0-foss-2018a

bedtools genomecov -ibam ${genome_prefix}.sorted.bam -bga > ${genome_prefix}.bedgraph
cat ${genome_prefix}.bedgraph | sort -k1,1 -k2,2n > ${genome_prefix}.sorted.bedgraph

module load ucsc/359

bedGraphToBigWig ${genome_prefix}.sorted.bedgraph /scratch/jl03308/Ab10_Assembly/analysis/Pipeline/pseudomolecules/B73Ab10.pseudomolecules-v2.chrs.sizes ${genome_prefix}.sorted.bw
#Unaligned positions are defined as regions with fewer than 3 reads aligned, bigger than 10kb, and distance bigger than 1kb with another unaligned region
cat ${genome_prefix}.bedgraph |awk '{if($4<3){print$0}}' |cut -f1,2,3 | sort -k1,1 -k2,2n |bedtools merge -i - -d 1000 |awk '{if($3-$2>10000){print$0}}' > ${genome_prefix}.readgap.bed

