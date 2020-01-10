#!/bin/bash

chip="$1"
control="$2"
genome="$3"
genome_name=$(basename $genome |cut -f1 -d ".")

module load SAMtools/1.9-foss-2016b

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
