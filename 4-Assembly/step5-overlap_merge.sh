'''Overlap merging step'''
module load Anaconda3/2019.03
module load Biopython/1.71-foss-2018a-Python-3.6.4 

xmap=$1
key=$2
out_dir=$3
python Overlap_merge.py -x ./alignref/core_sup_merged_hybridassembly0_CTTAAG_0kb_0labels.xmap \
-k ./core_sup_merged_hybridassembly0_CTTAAG_0kb_0labels_key.txt \
-o ./overlap_merge_src

singularity exec /usr/local/singularity-images/seqkit-0.10.1.simg seqkit grep -v -f modified_seq.txt core_sup_merged_hybridassembly0.fasta > core_sup_merged_hybridassembly0.tmp.fasta
tmp=$(head -1 ${output_dir}/output_seq.txt)
cat core_sup_merged_hybridassembly0.tmp.fasta $tmp > core_sup_merged_hybridassembly0.rmoverlap.fasta

