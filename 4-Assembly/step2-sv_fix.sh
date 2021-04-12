#!/bin/bash
if [ "$#" -ne 3 ] ; then
echo "Files required:"
echo -e "\t\t(1) sequence contigs in fasta format (ngs contigs for SV calling)"
echo -e "\t\t(2) bionano de novo assembly in cmap format"
echo -e "\t\t(3) output direcotry"
echo "step2-sv_fix.sh <contigs.fasta> <bionano.cmap> <output_dir>" ;
exit 0
fi

module load BionanoSolve/3.4-06042019-foss-2016b

genome="$1"
bionano_genome="$2"
output_dir="$3"

genome_name=$(basename ${genome%.*})
genome_dirc=$(dirname ${genome})
genome_cmap=$(basename ${genome%.*})_CTTAAG_0kb_0labels.cmap
errbin=$(basename ${genome%.*})_CTTAAG_0kb_0labels.errbin
#FASTA to CMAP conversion (1 minute and 18 seconds)
perl $EBROOTBIONANOSOLVE/HybridScaffold/06042019/scripts/fa2cmap_multi_color.pl -i ${genome} -e CTTAAG 1

#SV calling, use adjusted cmap as reference, to identify indels in contigs after conflict resolving

cd ${output_dir}
python $EBROOTBIONANOSOLVE/Pipeline/06042019/runCharacterize.py \
-t $EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/RefAligner \
-r $bionano_genome \
-q ${genome_dirc}/${genome_cmap} \
-p $EBROOTBIONANOSOLVE/Pipeline/06042019 \
-a $EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/optArguments_nonhaplotype_noES_DLE1_saphyr.xml \
-n 12 > ${genome_name}.stats

python $EBROOTBIONANOSOLVE/Pipeline/06042019/runSV.py \
-t $EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/RefAligner \
-r $bionano_genome \
-q ${genome_dirc}/${genome_cmap} \
-o ${output_dir}/SV_identified \
-E ${output_dir}/alignref/$errbin \
-a $EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/optArguments_nonhaplotype_noES_DLE1_saphyr.xml \
-j 12 \
-T 12 

