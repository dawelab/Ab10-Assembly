#!/bin/bash
reads=$1
out_dir=$2

ml minimap2
ml miniasm/0.3-foss-2016b 
ml Racon

minimap2 -x ava-ont -t 64 $reads $reads | gzip -1 > $out_dir/ont.allvsall.paf.gz
miniasm -f $reads $out_dir/ont.allvsall.paf.gz > $out_dir/ont.allvsall.gfa
awk '/^S/{print ">"$2"\n"$3}' $out_dir/ont.allvsall.gfa |fold > $out_dir/ont.rawcontigs.fa

minimap $out_dir/ont.rawcontigs.fa $reads > $out_dir/ont.ccs.1.paf
racon -t 64 $reads $out_dir/ont.ccs.1.paf $out_dir/ont.rawcontigs.fa $out_dir/ont.ccs.1.fa

minimap $out_dir/ont.ccs.1.fa $reads > $out_dir/ont.ccs.2.paf
racon -t 64 $reads $out_dir/ont.ccs.2.paf $out_dir/ont.ccs.1.fa $out/ont.ccs.2.fa

minimap $out_dir/ont.ccs.2.fa $reads > $out_dir/ont.ccs.3.paf
racon -t 64 $reads $out_dir/ont.ccs.3.paf $out_dir/ont.ccs.2.fa $out/ont.ccs.3.fa

