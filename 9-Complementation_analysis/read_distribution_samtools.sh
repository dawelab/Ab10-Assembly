#!/bin/bash

bam_file="$1"
region="$2"

module load SAMtools/1.9-foss-2016b
samtools view $bam_file $region |gawk '{print length($10)}' |awk '{if($1>10000){print$0}}' > readdistribution.txt


