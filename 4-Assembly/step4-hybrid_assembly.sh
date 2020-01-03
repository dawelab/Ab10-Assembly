#!/bin/bash

module load BionanoSolve/3.4-06042019-foss-2016b

ngs_genome=$1
bionano=$2
errbin=$3
outdir=$4
mkdir -p ${outdir}

perl $EBROOTBIONANOSOLVE/HybridScaffold/06042019/hybridScaffold.pl \
-n ${ngs_genome} \
-b ${bionano} \
-c $EBROOTBIONANOSOLVE/HybridScaffold/06042019/hybridScaffold_DLE1_config.xml \
-r $EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/sse/RefAligner \
-o ${outdir} \
-f \
-B 2 \
-N 2 \
-y \
-x \
-p $EBROOTBIONANOSOLVE/Pipeline/06042019 \
-q $EBROOTBIONANOSOLVE/RefAligner/8949.9232rel/optArguments_nonhaplotype_noES_noCut_DLE1_saphyr.xml \
-e ${errbin}
