#!/bin/bash
#PBS -q highmem_q
#PBS -N ${genome_name}.Input_centromere
#PBS -l nodes=1:ppn=12
#PBS -l walltime=25:00:00
#PBS -l mem=120gb
#PBS -M jl03308@uga.edu
#PBS -m ae

input_R1 = "$1"
input_R2 = "$2"
genome="$3"

genome_name=$(basename $genome |cut -f1 -d ".")

module load Trim_Galore/0.4.5-foss-2016b
trim_galore --fastqc --gzip --paired ${input_R1} ${input_R2} -o /scratch/jl03308/NAM_pancentromere/trimmed/Input/${genome_name}

module load BWA/0.7.17-foss-2016b
mkdir /scratch/jl03308/NAM_pancentromere/analysis/peak_call/${genome_name}
cd /scratch/jl03308/NAM_pancentromere/analysis/peak_call/${genome_name}
index=${genome%.%*}.ann
if [ -f "${index}" ]; then
   echo "skip genome indexing."
else
  bwa index $genome
fi

bwa mem -t 12 ${genome} /scratch/jl03308/NAM_pancentromere/trimmed/Input/${genome_name}/*R1_val_1.fq.gz /scratch/jl03308/NAM_pancentromere/trimmed/Input/${genome_name}/*_R2_val_2.fq.gz > ${genome_name}.Input.sam

module load SAMtools/1.9-foss-2016b
samtools view -@ 12 -b -o ${genome_name}.Input.bam ${genome_name}.Input.sam
samtools sort -o ${genome_name}.Input.sorted.bam -T ${genome_name}.Input -@ 12 ${genome_name}.Input.bam

ml picard/2.16.0-Java-1.8.0_144
java -jar /usr/local/apps/eb/picard/2.16.0-Java-1.8.0_144/picard.jar MarkDuplicates I=${genome_name}.Input.sorted.bam O=${genome_name}.Input.rmdup.bam M=${genome_name}.Input.marked_dup_metrics.txt REMOVE_DUPLICATES=true
samtools sort -o ${genome_name}.Input.rmdup.sorted.bam -T ${genome_name}.Input -@ 12 ${genome_name}.Input.rmdup.bam
samtools index -@ 12 ${genome_name}.Input.rmdup.sorted.bam
samtools view -@ 12 -bhq 20 ${genome_name}.Input.rmdup.sorted.bam -o ${genome_name}.Input.rmdup.q20.bam
samtools sort -o ${genome_name}.Input.rmdup.q20.sorted.bam -T ${genome_name}.Input -@ 12 ${genome_name}.Input.rmdup.q20.bam
samtools flagstat ${genome_name}.Input.rmdup.q20.sorted.bam > ${genome_name}.Input.flagstat
samtools index -@ 12 ${genome_name}.Input.rmdup.q20.sorted.bam

