#!/bin/bash
echo "Files required:"
echo -e "\t\t(1) contigs in fasta format (core/backbone assembly)"
echo -e "\t\t(2) contigs in fasta format (supplementary assembly)"
echo -e "\t\t(3) bionano contigs in cmap format (bionano de novo assembly)"
echo -e "\t\t(4) output directory for merged contigs and bionano alignment"
echo "";
echo "step3-merging.sh <core.fasta> <sup.cmap> <bionano.cmap> <output_dir>" ;

ml minimap2/2.13-foss-2016b 
module load miniasm/0.3-foss-2016b 
module load SAMtools/1.9-foss-2016b

core_genome="$1"
sup_genome="$2"
bionano_genome="$3"
output_dir="$4"

'''Merge core and supplementary genome through Miniasm'''

cd $output_dir
cat $core_genome $sup_genome > core_sup_all.fasta

minimap2 -t22 -k28 -w28 -A1 -B9 -O16,41 -E2,1 -z200 -g100000 -r100000 --max-chain-skip 100 $core_genome $sup_genome > core_sup_merged.paf
miniasm -1 -2 -r0 -e1 -n1 -I0.8 -h250000 -g100000 -o25000 core_sup_merged.paf > core_sup_merged_noseq.gfa
miniasm -1 -2 -r0 -e1 -n1 -I0.8 -h250000 -g100000 -o25000 core_sup_merged.paf -f core_sup_all.fasta > core_sup_merged.gfa
awk '/^S/{print ">"$2"\n"$3}' core_sup_merged.gfa | fold > core_sup_merged.fasta

# check nucleotides incorporated in the merged assembly
cat core_sup_merged_noseq.gfa | awk '{if($1=="a"&&$4~"core"){print$0}}' |cut -f6 | awk '{ SUM += $1} END { print SUM}' # core 
cat core_sup_merged_noseq.gfa | awk '{if($1=="a"&&$4!~"core"){print$0}}' |cut -f6 | awk '{ SUM += $1} END { print SUM}' # sup 
# Extract core contigs unincorporated in the merged assembly
cat core_sup_merged_noseq.gfa | awk '{if($1=="a"&&$4~"core"){print$0}}' |cut -f4 > core_assembled.ID 
singularity exec /usr/local/singularity-images/seqkit-0.10.1.simg seqkit grep -v -f core_assembled.ID $core_genome > core_SV_fixed.unassembled.fasta
cat core_sup_merged.fasta core_SV_fixed.unassembled.fasta > core_sup_merged_hybridassembly0.fasta # used for hybrid assembly
samtools faidx core_sup_merged_hybridassembly0.fasta 

'''Overlap merging step0'''
# alignment to bionano maps
perl $BionanoSolve/HybridScaffold/06042019/scripts/fa2cmap_multi_color.pl -i core_sup_merged_hybridassembly0.fasta -e CTTAAG 1 

python $BionanoSolve/Pipeline/06042019/runCharacterize.py \
-t $BionanoSolve/RefAligner/8949.9232rel/RefAligner \
-r $bionano_genome \
-q core_sup_merged_hybridassembly0_CTTAAG_0kb_0labels.cmap \
-p $BionanoSolve/Pipeline/06042019 \
-a $BionanoSolve/RefAligner/8949.9232rel/optArguments_nonhaplotype_noES_DLE1_saphyr.xml \
-n 12 

