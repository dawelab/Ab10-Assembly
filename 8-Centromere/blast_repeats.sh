#!/bin/bash

genome="$1"
genome_name=$(basename $genome |cut -f1 -d ".")
genome_dir=$(dirname $genome)
out_dir="$2"

cd ${out_dir}

module load BLAST/2.2.26-Linux_x86_64
formatdb -p F -i $genome
blastall -p blastn -d $genome -i ~/reference/knob180.fa -m 8 -b 5000 > knob180.blast
blastall -p blastn -d $genome -i ~/reference/TR-1.fa -m 8 -b 5000 > TR-1.blast
blastall -p blastn -d $genome -i ~/reference/CentC.fa -m 8 -b 5000 > CentC.blast
blastall -p blastn -d $genome -i ~/reference/rDNA_intergenic_spacer.fa -m 8 -b 5000 > rDNA_intergenic_spacer.blast
blastall -p blastn -d $genome -i ~/reference/subtelomeric_4-12-1.fa -m 8 -b 5000 > subtelomeric_4-12-1.blast
blastall -p blastn -d $genome -i ~/reference/Kindr.fa -m 8 -b 5000 > Kindr.blast
cat *.blast > blast.sum

module load seqtk 
seqtk cutN -n 10 -g $genome > ${genome_dir}/${genome_name}.Ngap.bed

module load Anaconda3/2019.03
module load Biopython/1.71-foss-2018a-Python-3.6.4 
python ~/repeat_analyses.py --blast_file ${out_dir}/blast.sum --gap_file ${genome_dir}/${genome_name}.Ngap.bed --output_dir {out_dir}
