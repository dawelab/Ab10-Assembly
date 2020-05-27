#!/bin/bash

TE_annotation="$1" 
nonunique="$2"

dir=$(dirname $nonunique)
cd $dir

line=$(basename $TE_annotation|cut -f1 -d ".")
#Extract transposon positions in genome
cat $TE_annotation |awk '{if($9~"huck"){print$0}}' |cut -f1,4,5 |cut -f2,3 -d "_" |sort -k1,1 -k2,2n  > ${line}.huck.bed
cat $TE_annotation |awk '{if($9~"cinful_zeon"){print$0}}' |cut -f1,4,5 |cut -f2,3 -d "_" |sort -k1,1 -k2,2n  > ${line}.cinful-zeon.bed
cat $TE_annotation |awk '{if($9~"grande"){print$0}}' |cut -f1,4,5 |cut -f2,3 -d "_" |sort -k1,1 -k2,2n  > ${line}.grande.bed
cat $TE_annotation |awk '{if($9~"prem1"||$9~"xilon_diguus"||$9~"tekay"){print$0}}' |cut -f1,4,5 |cut -f2,3 -d "_" |sort -k1,1 -k2,2n  > ${line}.prem1.bed
cat $TE_annotation |awk '{if($9~"opie"||$9~"ji"||$9~"ruda"||$9~"giepum"){print$0}}' |cut -f1,4,5 |cut -f2,3 -d "_" |sort -k1,1 -k2,2n  > ${line}.opie-ji.bed
cat $TE_annotation |awk '{if($9~"CRM"){print$0}}' |cut -f1,4,5 |cut -f2,3 -d "_" |sort -k1,1 -k2,2n  > ${line}.crm.bed
cat $TE_annotation |awk '{if($9~"CentC"){print$0}}' |cut -f1,4,5 |cut -f2,3 -d "_" |sort -k1,1 -k2,2n  > ${line}.centC.bed


module load BEDTools/2.28.0-foss-2018a

#step1:find TEs intersect with active centromere; step2: divide the TEs into unique/nonunique groups (based on intersection status with nonunique regions)
for file in ${line}.*.bed
do
type=$(basename $file | cut -f2 -d ".")
# the intersection of transposons with active centromeres reports all transposons that have overlaps with centromere coords, not necessarily part of a transposon on the boundary of a centromere 
bedtools intersect -u -nonamecheck -a $file -b ${line}.active_centromere | sort -k1,1 -k2,2n |bedtools intersect -f 0.1 -u -nonamecheck -a - -b $nonunique |sort -k1,1 -k2,2n > ${line}.nonunique_centromere.${type}.bed 
bedtools intersect -u -nonamecheck -a $file -b ${line}.active_centromere | sort -k1,1 -k2,2n |bedtools intersect -f 0.1 -v -nonamecheck -a - -b $nonunique |sort -k1,1 -k2,2n > ${line}.unique_centromere.${type}.bed 
done

for file in ${line}.*unique_centromere.*.bed
do
type=$(basename $file |cut -f1,2,3 -d ".")
bedtools map -a $file -b ${line}.ChIP_Input.RPKM.sorted.bedgraph -c 4 -o mean > ${type}.RPKM.bedgraph
done

#Quantify the abundance and enrichment of uniquely and nonuniquely mapped centc and 6 transposon families.
for type in $(ls ${line}.nonunique_centromere.*.bed  |cut -f3 -d "." |sort |uniq); do 
for x in chr{1..10}; do

TE_abundance=$(bedtools intersect -u -nonamecheck -b ${line}.active_centromere -a ${line}.${type}.bed |sort -k1,1 -k2,2n |awk -v chr=$x '{if($1==chr){print$0}}'|awk '{ sum += $3-$2} END {print sum}')
unique_abundance=$(cat ${line}.unique_centromere.${type}.bed | awk -v chr=$x '{if($1==chr){print$0}}'| awk '{ sum += $3-$2} END {print sum}')
nonunique_abundance=$(cat ${line}.nonunique_centromere.${type}.bed | awk -v chr=$x '{if($1==chr){print$0}}'| awk '{ sum += $3-$2} END {print sum}')

unique_total_enrichement=$(cat ${line}.unique_centromere.${type}.RPKM.bedgraph  |awk -v chr=$x '{if($1==chr){print$0}}' |awk '{print$0"\t"$4*($3-$2)}'|awk '{ sum += $5 } END {print sum}') 
nonunique_total_enrichement=$(cat ${line}.nonunique_centromere.${type}.RPKM.bedgraph |awk -v chr=$x '{if($1==chr){print$0}}' |awk '{print$0"\t"$4*($3-$2)}'|awk '{ sum += $5 } END {print sum}') 

if [ -z "$unique_abundance" ]; then
unique_enrichement="NA"
else
unique_enrichement=$(awk -v a=$unique_total_enrichement -v b=$unique_abundance 'BEGIN{print (a / b)}')
fi

if [ -z "$nonunique_abundance" ]; then
nonunique_enrichement="NA"
else
nonunique_enrichement=$(awk -v a=$nonunique_total_enrichement -v b=$nonunique_abundance 'BEGIN{print (a / b)}')
fi

echo -e $line"\t"$type"\t"$x"\t"$TE_abundance"\t"$unique_abundance"\t"$unique_enrichement"\t"$nonunique_abundance"\t"$nonunique_enrichement"\t" >> ${line}.enrichment.sum.txt
done
done
