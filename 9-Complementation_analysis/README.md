We studied how nanopore and pacbio assemblies complemented each other in repetitive and heterozygous area. The complementation examination was performed by comparing the relative position of ont and pacbio contig assembly gaps to the final assembly. Length distribution of nanopore and pacbio reads mapped to these gaps were compared. 

Tandemly repeated areas by chromosome self-alignment with minimap2 (v2.17; -PD -k19 -w19 -m200) and heterozygous regions by manual inspection using Bionano Access software.

PacBio gap coordinates were projected onto the final assembly using minimap2 (v2.17; -cx asm5 --cs), followed by coordinate liftover using paftools.js. Gaps that were complemented by Nanopore contigs were identified as gaps present in the PacBio assembly but absent in the final assembly. The PacBio adjusted gap coordinates, complemented gaps, and final assembly gaps were mapped to tandem repeats and heterozygous regions with bedtools (v2.28.0; window -r 500000 -l 500000). The co-occurence of PacBio gaps with tandem repetitiveness and heterozygous regions was assessed by two-tailed Fisher’s exact test using bedtools fisher (v2.28.0) at default settings.


