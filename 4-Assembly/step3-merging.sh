#!/bin/bash
if [ "$#" -ne 4 ] ; then
echo "Files required:"
echo -e "\t\t(1) contigs in fasta format (core/backbone assembly)"
echo -e "\t\t(2) contigs in fasta format (supplementary assembly)"
echo -e "\t\t(3) bionano contigs in cmap format (bionano de novo assembly)"
echo -e "\t\t(4) output directory for merged contigs and bionano alignment"
echo "step3-merging.sh <core.fasta> <sup.fasta> <bionano.cmap> <output_dir>" ;
exit 0
fi

ml minimap2/2.13-foss-2016b 
module load miniasm/0.3-foss-2016b 
module load SAMtools/1.9-foss-2016b

core_genome="$1"
sup_genome="$2"
bionano_genome="$3"
output_dir="$4"

sed 's/>.*/&_core/' $core_genome > ${core_genome%.*}.renamed.fasta

'''Merge core and supplementary genome through Miniasm'''

cat ${core_genome%.*}.renamed.fasta $sup_genome > ${output_dir}/core_sup_all.fasta

minimap2 -t22 -k28 -w28 -A1 -B9 -O16,41 -E2,1 -z200 -g100000 -r100000 --max-chain-skip 100 ${core_genome%.*}.renamed.fasta $sup_genome > ${output_dir}/core_sup_merged.paf
miniasm -1 -2 -r0 -e1 -n1 -I0.8 -h250000 -g100000 -o25000 ${output_dir}/core_sup_merged.paf > ${output_dir}/core_sup_merged_noseq.gfa
miniasm -1 -2 -r0 -e1 -n1 -I0.8 -h250000 -g100000 -o25000 ${output_dir}/core_sup_merged.paf -f ${output_dir}/core_sup_all.fasta > ${output_dir}/core_sup_merged.gfa
awk '/^S/{print ">"$2"\n"$3}' ${output_dir}/core_sup_merged.gfa | fold > ${output_dir}/core_sup_merged.fasta

# check nucleotides incorporated in the merged assembly
cat ${output_dir}/core_sup_merged_noseq.gfa | awk '{if($1=="a"&&$4~"core"){print$0}}' |cut -f6 | awk '{ SUM += $1} END { print SUM}' # core 
cat ${output_dir}/core_sup_merged_noseq.gfa | awk '{if($1=="a"&&$4!~"core"){print$0}}' |cut -f6 | awk '{ SUM += $1} END { print SUM}' # sup 
# Extract core contigs unincorporated in the merged assembly
cat ${output_dir}/core_sup_merged_noseq.gfa | awk '{if($1=="a"&&$4~"core"){print$0}}' |cut -f4 > ${output_dir}/core_assembled.ID 
singularity exec /usr/local/singularity-images/seqkit-0.10.1.simg seqkit grep -v -f ${output_dir}/core_assembled.ID ${core_genome%.*}.renamed.fasta > ${output_dir}/core_SV_fixed.unassembled.fasta
cat ${output_dir}/core_sup_merged.fasta ${output_dir}/core_SV_fixed.unassembled.fasta > ${output_dir}/core_sup_merged_hybridassembly0.fasta # used for hybrid assembly
samtools faidx ${output_dir}/core_sup_merged_hybridassembly0.fasta 

'''Overlap merging step0'''

module load BionanoSolve/3.4-06042019-foss-2016b
# alignment to bionano maps
perl $EBROOTBIONANOSOLVE/HybridScaffold/06042019/scripts/fa2cmap_multi_color.pl -i ${output_dir}/core_sup_merged_hybridassembly0.fasta -e CTTAAG 1 -o ${output_dir}

python $EBROOTBIONANOSOLVE/Pipeline/06042019/runCharacterize.py \
-t $EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/RefAligner \
-r $bionano_genome \
-q ${output_dir}/core_sup_merged_hybridassembly0_CTTAAG_0kb_0labels.cmap \
-p $EBROOTBIONANOSOLVE/Pipeline/06042019 \
-a $EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/optArguments_nonhaplotype_noES_DLE1_saphyr.xml \
-n 12 
-o ${output_dir}

