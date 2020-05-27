#!/bin/bash

bam_input="$1" 
bam_chip="$2" 

bam_dir=$(dirname $bam_input)
line=$(basename $bam_input|cut -f1 -d ".")

cd $bam_dir

#Identify the enrichment (RPKM ratio/log2 ratio) in genome
module load Python/3.6.4-foss-2018a
module load deepTools/3.2.1-foss-2018a-Python-3.6.4

bamCompare --bamfile1 $bam_chip --bamfile2 $bam_input --binSize 1 --scaleFactorsMethod None --normalizeUsing RPKM --numberOfProcessors max --ignoreDuplicates --minMappingQuality 20 --outFileFormat bedgraph --outFileName ${line}.ChIP_Input.RPKM.log.bedgraph #log ratio
bamCompare --bamfile1 $bam_chip --bamfile2 $bam_input --binSize 1 --operation ratio --scaleFactorsMethod None --normalizeUsing RPKM --numberOfProcessors max --ignoreDuplicates --minMappingQuality 20 --outFileFormat bedgraph --outFileName ${line}.ChIP_Input.RPKM.bedgraph #ratio
cat ${line}.ChIP_Input.RPKM.bedgraph |sort -k1,1 -k2,2n > ${line}.ChIP_Input.RPKM.sorted.bedgraph #sort bedgraph file
