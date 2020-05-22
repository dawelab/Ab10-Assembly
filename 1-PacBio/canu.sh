#!/bin/bash
output_dir=$1
reads=$1

ml canu/1.8-Linux-amd64
canu -d $output -p ab10 genomeSize=2.3g -pacbio-corrected $reads gridOptions="-V " correctedErrorRate=0.065 corMhapSensitivity=normal ovlMerThreshold=500 utgOvlMerThreshold=150
