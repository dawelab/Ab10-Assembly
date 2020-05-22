#!/bin/bash -login
#SBATCH -N 1
#SBATCH -t 7-0:00:00
#SBATCH --ntasks-per-node 36
##SBATCH --mem=300GB

echo 'Start:'
date

conda activate EDTA

mkdir B73Ab10
cd B73Ab10

#download the curated MTEC library
wget https://raw.githubusercontent.com/oushujun/MTEC/master/maizeTE02052020

#annotate the genome with EDTA
/usr/bin/time -v perl ~/las/git_bin/EDTA/EDTA.pl --genome B73Ab10.pseudomolecules-v2.fasta --species Maize -t 36 --anno 1 --curatedlib maizeTE02052020 --cds zea_maysb73ab10.evd.cds.fasta

echo 'End:'
date

