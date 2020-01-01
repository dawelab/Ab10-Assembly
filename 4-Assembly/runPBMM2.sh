#!/bin/bash
if [ "$#" -ne 2 ] ; then
echo "Please provide a genome file and the subread (bam) file";
echo "./runPBMM2.sh <genome.fasta> <subreads.bam>"
exit 0
fi
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate pacbio
ref=$1
subreads=$2
index=${ref%.*}.mmi
if [ -f "${index}" ]; then
   echo "mmi file found, skipping genome fasta indexing.."
else
  pbmm2 index ${ref} ${index}
fi
out=$(basename ${subreads%.*})
pbmm2 align ${ref} ${subreads} ${out}-Ab10.bam
ml samtools
samtools sort -o ${out}-Ab10_sorted.bam -T ${out}.temp --threads 36 ${out}-Ab10.bam
