module load Anaconda3/2019.03
module load Biopython/1.71-foss-2018a-Python-3.6.4 

xmap=$1
key=$2
out_dir=$3
python Overlap_merge.py -x $xmap \
-k $key \
-o $out_dir

