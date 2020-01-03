#!/bin/bash

module load BionanoSolve/3.4-06042019-foss-2016b

genome="$1"
bionano_genome="$2"
output_dir="$3"

bionano_name=$(basename ${bionano_genome%.*})
genome_name=$(basename ${genome%.*})
genome_dirc=$(dirname ${genome})
genome_cmap=$(basename ${genome%.*})_CTTAAG_0kb_0labels.cmap
errbin=$(basename ${genome%.*})_CTTAAG_0kb_0labels.errbin
#FASTA to CMAP conversion (1 minute and 18 seconds)
perl $EBROOTBIONANOSOLVE/HybridScaffold/06042019/scripts/fa2cmap_multi_color.pl -i ${genome} -e CTTAAG 1

# Cut conflicts for core assembly 
'''Conflict-resolving on cmap level'''
# Align0, generate adjusted cmap 
$EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/RefAligner \
-ref ${genome_dirc}/${genome_cmap} \
-i $bionano_genome \
-o ${output_dir}/align0 \
-stdout -stderr -maxmem 128 -maxthreads 22 -RAmem 3 1 -M 1 3 -ScaleDelta 0.02 15 \
-ScaleDeltaBPP -hashScaleDelta 2 -res 2.9 -FP 0.6 -FN 0.06 -sf 0.20 \
-sd 0.0 -sr 0.01 -extend 1 -outlier 0.0001 -endoutlier 0.001 -PVendoutlier -deltaX 12 -deltaY 12 \
-xmapchim 12 -hashgen 5 7 2.4 1.5 0.05 5.0 1 1 4 -hash -hashdelta 26 10 46 -hashMultiMatch 30 10 -insertThreads 4 -nosplit 2 -biaswt 0 -T 1e-10 -S -1000 \
-indel -PVres 2 -rres 0.9 -MaxSE 0.5 -HSDrange 1.0 -outlierBC -xmapUnique 12 -AlignRes 2. -outlierExtend 12 24 -Kmax 12 -resEstimate -f -mres 0.9

$EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/RefAligner \
-merge -i $bionano_genome \
-o ${output_dir}/${bionano_name}_bppAdjust \
-readparameters ${output_dir}/align0.errbin -stdout -stderr

# Align1

$EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/RefAligner \
-ref ${genome_dirc}/${genome_cmap} \
-i ${output_dir}/${bionano_name}_bppAdjust.cmap \
-o ${output_dir}/align1 -stdout -stderr -maxmem 128 -maxthreads 64 -maxvirtmem 0 -RAmem 3 1 -res 2.9 -FP 0.6 -FN 0.06 -sf 0.20 -sd 0.0 -sr 0.01 -extend 1 \
-outlier 0.0001 -endoutlier 0.001 -PVendoutlier -deltaX 12 -deltaY 12 -xmapchim 12 -hashgen 5 7 2.4 1.5 0.05 5.0 1 1 4 -hash -hashdelta 26 10 46 \
-hashMultiMatch 30 10 -insertThreads 4 -nosplit 2 -biaswt 0 -T 1e-10 -S -1000 -indel -PVres 2 -rres 0.9 -MaxSE 0.5 -HSDrange 1.0 \
-outlierBC -xmapUnique 12 -AlignRes 2. -outlierExtend 12 24 -Kmax 12 -resEstimate -f -mres 0.9

# AlignAssignType

perl $EBROOTBIONANOSOLVE/HybridScaffold/06042019/scripts/AssignAlignType.pl \
${output_dir}/align1.xmap \
${output_dir}/align1_r.cmap \
${output_dir}/align1_q.cmap \
${output_dir}/assignAlignType.xmap \
${output_dir}/assignAlignType_r.cmap \
${output_dir}/assignAlignType_q.cmap 11 5 \
${genome_dirc}/${genome_cmap} \
${output_dir}/${bionano_name}_bppAdjust.cmap \
${output_dir}/conflicts.txt 30

# Cut conflicts

perl $EBROOTBIONANOSOLVE/HybridScaffold/06042019/scripts/cut_conflicts.pl \
-align1XmapFile ${output_dir}/align1.xmap \
-align1GMFile ${output_dir}/align1_q.cmap \
-align1SeqFile ${output_dir}/align1_r.cmap \
-maxOverhang 10 -breakPointFileShiftAmount 30 \
-oriGMFile ${output_dir}/${bionano_name}_bppAdjust.cmap \
-oriSeqFile ${genome_dirc}/${genome_cmap} \
-conflictFile ${output_dir}/conflicts.txt \
-outDir ${output_dir}/conflict_CMAP \
-outGMFile ${output_dir}/conflict_CMAP/${bionano_name}_bppAdjust_cut.cmap \
-outSeqFile ${output_dir}/conflict_CMAP/$(basename ${genome%.*})_CTTAAG_0kb_0labels_cut.cmap \
-windowSize 10000 -qScoreThreshold 35 -covThreshold 10 \
-refAligner $EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/RefAligner

'''Conflict-resolving on sequence level'''
module load Anaconda3/5.0.1
module load Biopython/1.68-foss-2016b-Python-3.5.2

mkdir ${output_dir}/conflict_NGS
python ~/cut_conflicts.py \
-i $genome \
-k $core/step0_assembly/${genome_name}_CTTAAG_0kb_0labels_key.txt \
-c ${output_dir}/conflict_CMAP/auto_cut_NGS_coord_translation.txt \
-o ${output_dir}/conflict_NGS/${genome_name}_cut.fasta

