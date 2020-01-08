#!/bin/bash

chip_SE="$1"
genome="$2"
genome_name=$(basename $genome |cut -f1 -d ".")

module load Trim_Galore/0.4.5-foss-2016b
trim_galore --fastqc --gzip $chip_SE -o /scratch/jl03308/NAM_pancentromere/trimmed/ChIP/${genome_name}

module load BWA/0.7.17-foss-2016b
cd /scratch/jl03308/NAM_pancentromere/analysis/peak_call/${genome_name}
index=${genome%.%*}.ann
if [ -f "${index}" ]; then
   echo "skip genome indexing."
else
  bwa index $genome
fi
bwa mem -t 12 ${genome} /scratch/jl03308/NAM_pancentromere/trimmed/ChIP/${genome_name}/*.fq.gz > ${genome_name}.ChIP.sam

module load SAMtools/1.9-foss-2016b
samtools view -@ 12 -b -o ${genome_name}.ChIP.bam ${genome_name}.ChIP.sam
samtools sort -o ${genome_name}.ChIP.sorted.bam -T ${genome_name}.ChIP -@ 12 ${genome_name}.ChIP.bam

ml picard/2.16.0-Java-1.8.0_144
java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates I=${genome_name}.ChIP.sorted.bam O=${genome_name}.ChIP.rmdup.bam M=${genome_name}.ChIP.marked_dup_metrics.txt REMOVE_DUPLICATES=true
samtools sort -o ${genome_name}.ChIP.rmdup.sorted.bam -T ${genome_name}.ChIP -@ 12 ${genome_name}.ChIP.rmdup.bam
samtools index -@ 12 ${genome_name}.ChIP.rmdup.sorted.bam
samtools view -@ 12 -bhq 20 ${genome_name}.ChIP.rmdup.sorted.bam -o ${genome_name}.ChIP.rmdup.q20.bam
samtools sort -o ${genome_name}.ChIP.rmdup.q20.sorted.bam -T ${genome_name}.ChIP -@ 12 ${genome_name}.ChIP.rmdup.q20.bam
samtools flagstat ${genome_name}.ChIP.rmdup.q20.sorted.bam > ${genome_name}.ChIP.flagstat
samtools index -@ 12 ${genome_name}.ChIP.rmdup.q20.sorted.bam

