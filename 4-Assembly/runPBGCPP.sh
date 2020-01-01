#!/bin/bash
if [ "$#" -ne 2 ] ; then
echo "please provide a genome (fasta) and bam file (subreads algined to genome)";
echo "./runPBGCPP.sh <genome.fasta> <subreads-mapped-to-genome.bam>
exit 0
fi
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate pacbio
genome=$1
bam=$2
gcpp -r $genome -o ${genome%.*}-gcpp-polished.fasta --log-file polish-log.txt ${bam}
