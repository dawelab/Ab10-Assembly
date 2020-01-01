#!/bin/bash
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate repeatmasker
lib="/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/b-repeatmasking/TE_12-Feb-2015_15-35.fa"
genome=$1
RepeatMasker \
  -e ncbi \
  -pa 36 \
  -q \
  -lib ${lib} \
  -nocut \
  -gff \
  ${genome}

