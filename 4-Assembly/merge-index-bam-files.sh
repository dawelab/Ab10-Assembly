#!/bin/bash
ls *_sorted.bam > bamfiles.fofn
ml samtools
samtools merge --threads 36 -O BAM -b bamfiles.fofn B73Ab10_subreads.bam
samtools sort -o B73Ab10_subreads_sorted.bam -T temp-ab10-bam --threads 36 B73Ab10_subreads.bam
samtools index -@ 36 B73Ab10_subreads_sorted.bam
