#!/bin/bash

module load BionanoSolve/3.4-06042019-foss-2016b
module load SAMtools/1.9-foss-2016b
ml minimap2/2.13-foss-2016b 
module load miniasm/0.3-foss-2016b 


core_genome="$1"
sup_genome="$2"
bionano_genome="$3"
bionano_errbin="$4"
output_dir="$5"

core_name=$(basename ${core_genome%.*})

sup_name=$(basename ${sup_genome%.*})

mkdir -p ${output_dir}/core_assembly
mkdir -p ${output_dir}/core_assembly/step0_assembly
mkdir -p ${output_dir}/core_assembly/step1_cut_conflict
mkdir -p ${output_dir}/core_assembly/step2_sv_fix
mkdir -p ${output_dir}/core_assembly/step3_merging
mkdir -p ${output_dir}/core_assembly/step3_merging/contig_merge
mkdir -p ${output_dir}/core_assembly/step3_merging/contig_overlap_merge
mkdir -p ${output_dir}/core_assembly/step4_hybridassembly 
mkdir -p ${output_dir}/core_assembly/step4_hybridassembly/hybridassembly0 
mkdir -p ${output_dir}/core_assembly/step4_hybridassembly/final_scaffolds
mkdir -p ${output_dir}/core_assembly/step4_hybridassembly/overlap_merge

mkdir -p ${output_dir}/sup_assembly
mkdir -p ${output_dir}/sup_assembly/step0_assembly
mkdir -p ${output_dir}/sup_assembly/step1_cut_conflict
mkdir -p ${output_dir}/sup_assembly/step2_sv_fix
mkdir -p ${output_dir}/bionano_assembly

cp $core_genome ${output_dir}/core_assembly/step0_assembly/ .
core_genome_file=${output_dir}/core_assembly/step0_assembly/$(basename ${core_genome})
cp $sup_genome ${output_dir}/sup_assembly/step0_assembly/ .
sup_genome_file=${output_dir}/sup_assembly/step0_assembly/$(basename ${sup_genome})
cp $bionano_genome ${output_dir}/bionano_assembly/ .
cp $bionano_errbin ${output_dir}/bionano_assembly/ .

core_cut_output_dir="$5"
sup_cut_output_dir="$6"
core_sv_output_dir="$7"
sup_sv_output_dir="$8"

sh step1-conflict_resolving.sh $core_genome_file $bionano_genome ${output_dir}/core_assembly/step1_cut_conflict
sh step1-conflict_resolving.sh $sup_genome_file $bionano_genome ${output_dir}/sup_assembly/step1_cut_conflict

core_genome_cut=${output_dir}/core_assembly/step1_cut_conflict/conflict_NGS/${core_name}_cut.fasta
sup_genome_cut=${output_dir}/sup_assembly/step1_cut_conflict/conflict_NGS/${sup_name}_cut.fasta

sh step2-sv_fix.sh $core_genome_cut $bionano_genome ${output_dir}/core_assembly/step2_sv_fix
sh step2-sv_fix.sh $sup_genome_cut $bionano_genome ${output_dir}/sup_assembly/step2_sv_fix

module load Anaconda3/2019.03
module load Biopython/1.71-foss-2018a-Python-3.6.4 

python SV_fix.py \ 
-cs ${output_dir}/core_assembly/step2_sv_fix/SVs_identified \ 
-ss ${output_dir}/sup_assembly/step2_sv_fix/SVs_identified \
-cc ${output_dir}/core_assembly/step1_cut_conflict/conflict_CMAP/auto_cut_NGS_coord_translation.txt \
-sc ${output_dir}/sup_assembly/step1_cut_conflict/conflict_CMAP/auto_cut_NGS_coord_translation.txt \
-o ${output_dir}/core_assembly/step2_sv_fix

sh step3-merging.sh ${output_dir}/core_assembly/step2_sv_fix/core_SV_fixed.fasta $sup_genome_cut $bionano_genome ${output_dir}/core_assembly/step3_merging/contig_merge

#overlap_merge
python Overlap_merge.py \ 
-x ${output_dir}/core_assembly/step3_merging/contig_merge/alignref/core_sup_merged_hybridassembly0_CTTAAG_0kb_0labels.xmap \
-k ${output_dir}/core_assembly/step3_merging/contig_merge/core_sup_merged_hybridassembly0_CTTAAG_0kb_0labels_key.txt \
-o ${output_dir}/core_assembly/step3_merging/contig_overlap_merge

singularity exec /usr/local/singularity-images/seqkit-0.10.1.simg seqkit grep -v -f ${output_dir}/core_assembly/step3_merging/contig_overlap_merge/modified_seq.txt ${output_dir}/core_assembly/step3_merging/contig_merge/core_sup_merged_hybridassembly0.fasta > ${output_dir}/core_assembly/step3_merging/contig_merge/core_sup_merged_hybridassembly0.tmp.fasta
tmp=$(head -1 ${output_dir}/core_assembly/step3_merging/contig_overlap_merge/output_seq.txt)
cat ${output_dir}/core_assembly/step3_merging/contig_merge/core_sup_merged_hybridassembly0.tmp.fasta $tmp > ${output_dir}/core_assembly/step3_merging/contig_overlap_merge/core_sup_merged_hybridassembly0.rmoverlap.fasta

sh step4-hybrid_assembly.sh ${output_dir}/core_assembly/step3_merging/contig_overlap_merge/core_sup_merged_hybridassembly0.rmoverlap.fasta $bionano_genome $errbin ${output_dir}/core_assembly/step4_hybridassembly/hybridassembly0 

#overlap_merge_posthybridassembly
python Overlap_merge.py \ 
-x ${output_dir}/core_assembly/step4_hybridassembly/hybridassembly0/*_NGScontigs_HYBRID_SCAFFOLD.xmap \
-k ${output_dir}/core_assembly/step4_hybridassembly/hybridassembly0/*_key.txt.cut.txt \
-o ${output_dir}/core_assembly/step4_hybridassembly/overlap_merge

singularity exec /usr/local/singularity-images/seqkit-0.10.1.simg seqkit grep -v -f ${output_dir}/core_assembly/step4_hybridassembly/overlap_merge/modified_seq.txt ${output_dir}/core_assembly/step4_hybridassembly/hybridassembly0/*fa.cut.fasta > ${output_dir}/core_assembly/step4_hybridassembly/overlap_merge/core_sup_merged_hybridassembly1.tmp.fasta
tmp=$(head -1 ${output_dir}/core_assembly/step4_hybridassembly/overlap_merge/output_seq.txt)
cat ${output_dir}/core_assembly/step4_hybridassembly/overlap_merge/core_sup_merged_hybridassembly1.tmp.fasta $tmp > ${output_dir}/core_assembly/step4_hybridassembly/overlap_merge/core_sup_merged_hybridassembly1.rmoverlap.fasta

sh step4-hybrid_assembly.sh ${output_dir}/core_assembly/step4_hybridassembly/overlap_merge/core_sup_merged_hybridassembly1.rmoverlap.fasta $bionano_genome $errbin ${output_dir}/core_assembly/step4_hybridassembly/final_scaffolds 

