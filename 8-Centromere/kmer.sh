#!/bin/bash
#PBS -q highmem_q
#PBS -N kmer
#PBS -l nodes=1:ppn=12
#PBS -l walltime=1:30:00
#PBS -l mem=250gb
#PBS -M jl03308@uga.edu
#PBS -m ae

genome="$1"
genome_name=$(basename $genome |cut -f1 -d ".")
cd /scratch/jl03308/NAM_pancentromere/analysis/peak_call/${genome_name}

module load Jellyfish/2.2.6-foss-2016b

jellyfish count -t 12 -m 150 -s 2193M \
          --out-counter-len 1 --counter-len 1 $genome
jellyfish stats mer_counts.jf > ${genome_name}.stats
