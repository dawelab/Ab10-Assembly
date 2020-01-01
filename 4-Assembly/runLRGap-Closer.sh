#!/bin/bash
PATH=$PATH:/work/LAS/mhufford-lab/arnstrm/programs/LR_Gapcloser
ml perl-bio-perl
pacbio="/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/n_gapclosing/ec-subreads/pacbio-ec.fa"
ont="/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/n_gapclosing/ec-subreads/ab10.correctedReads.renamed.fasta"
#genome="/work/LAS/mhufford-lab/arnstrm/newNAM/analyses/n_gapclosing/all/EXP_REFINEFINAL1_bppAdjust_cmap_merged_2nd_fa_NGScontigs_HYBRID_SCAFFOLD.all.fasta"
genone="EXP_REFINEFINAL1_bppAdjust_cmap_merged_2nd_fa_NGScontigs_HYBRID_SCAFFOLD.all-gapclosed-pacbio.fasta"
#mkdir -p lrgapout-pcb
LR_Gapcloser.sh -i EXP_REFINEFINAL1_bppAdjust_cmap_merged_2nd_fa_NGScontigs_HYBRID_SCAFFOLD.all-gapclosed-pacbio.fasta -l /work/LAS/mhufford-lab/arnstrm/newNAM/analyses/n_gapclosing/ec-subreads/pacbio-ec.fa -s p -t 36 -o lrgapout-pcb2
