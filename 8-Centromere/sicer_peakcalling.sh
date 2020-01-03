#!/bin/bash
#PBS -q gpu_q
#PBS -N ${genome_name}.sicer
#PBS -l nodes=1:ppn=2
#PBS -l walltime=1:00:00
#PBS -l mem=30gb
#PBS -M jl03308@uga.edu
#PBS -m ae

chip="$1"
control="$2"
genome="$3"
genome_name=$(basename $genome |cut -f1 -d ".")

faidx=${genome%.%*}.fai
if [ -f "${faidx}" ]; then
   echo "skip genome faidx."
else
  samtools faidx ${genome}
fi

cat ${genome}.fai | cut -f1,2 > ${genome%.*}.chrom.sizes


module load epic2/0.0.36_conda
epic2 --treatment $chip --control $control --chromsizes ${genome%.*}.chrom.sizes --effective-genome-fraction 0.80 --bin-size 5000 --gaps-allowed 0 --output ${genome_name}.0gap_5kw.bed
cat ${genome_name}.0gap_5kw.bed|awk '{if($10>2&&$5>250){print$0}}' > ${genome_name}.0gap_5kw.filtered.bed