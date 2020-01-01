# Gene Prediction

## 1. Evidence Based Gene Prediction

Overview of steps:

The steps for evidence based prediction are as follows:

1. Map RNAseq reads to the genome
2. Assemble transcripts using various transcript assemblers
3. Run Mikado to consolidate transcripts and run programs to identify splice junctions, ORFs and full length matches to SwissProt plants proteins.
4. Pick transcripts for each locus and finalize evidence based annotations.
5. Post-processing of evidence based predictions.

Steps in detail:

### 1. Map RNAseq reads:

The RNAseq reads generated for Ab10 were as follows:

| ID      | PlantStage | Tissue      | Type   | TotalReads | AvgLength |
|---------|-------------|-------------|--------|------------:|-----------:|
| MN02011 | 8DAS        | root        | Paired |  29,818,892 |        152 |
| MN02012 | 8DAS        | root        | Paired |  32,543,108 |        152 |
| MN02021 | 8DAS        | shoot       | Paired |  26,420,582 |        152 |
| MN02022 | 8DAS        | shoot       | Paired |  32,247,423 |        152 |
| MN02031 | V11         | leaf-base   | Paired |  22,147,809 |        151 |
| MN02032 | V11         | leaf-base   | Paired |  27,560,950 |        151 |
| MN02041 | V11         | leaf-middle | Paired |  22,146,843 |        151 |
| MN02042 | V11         | leaf-middle | Paired |  22,466,821 |        151 |
| MN02051 | V11         | leaf-tip    | Paired |  26,407,717 |        151 |
| MN02052 | V11         | leaf-tip    | Paired |  21,923,046 |        151 |
| MN02061 | V18         | tassel      | Paired |  21,278,567 |        151 |
| MN02062 | V18         | tassel      | Paired |  22,868,224 |        151 |
| MN02071 | V18         | ear         | Paired |  30,154,390 |        151 |
| MN02072 | V18         | ear         | Paired |  20,563,411 |        151 |
| MN02081 | R1          | anther      | Paired |  17,818,298 |        151 |
| MN02082 | R1          | anther      | Paired |  24,882,545 |        151 |
| MN02091 | 16DAP       | endosperm   | Paired |  26,834,976 |        151 |
| MN02092 | 16DAP       | endosperm   | Paired |  26,182,541 |        151 |
| MN02101 | 16DAP       | embryo      | Paired |  14,833,469 |        151 |
| MN02102 | 16DAP       | embryo      | Paired |  26,947,220 |        151 |


Reads were mapped using STAR mapping program as follows:

```bash
# index
runSTARmap_index.sh B73Ab10.pseudomolecules-v2.fasta
# first round mapping
for read1 in *R1.fq.gz; do
  runSTARmap_round1.sh B73Ab10.pseudomolecules-v2.fasta ${read1}
done
# consolidate splice info
awk -f sjCollapseSamples.awk *_SJ.out.tab | sort -k1,1V -k2,2n -k3,3n > SJ.all
# second round mapping
for read1 in *R1.fq.gz; do
  runSTARmap_round2.sh B73Ab10.pseudomolecules-v2.fasta ${read1}
done
# merge BAM files
merge-bam-files.sh B73_AB10
```

### 2. Transcript Assembly

Five different transcript assemblers were ran on the merged BAM file (`merged_B73_AB10.bam`) as follows:


```bash
runStrawberry.sh merged_B73_AB10.bam
runStringtie.sh merged_B73_AB10.bam
runTrinity.sh merged_B73_AB10.bam
runClass2.sh merged_B73_AB10.bam
runCufflinks.sh merged_B73_AB10.bam
```

Since Trinity generates a fasta file, we will use GMAP to map them and create a GFF3 file.

```bash
runGMAPcdna.sh B73Ab10.pseudomolecules-v2.fasta
```

We will also need splice junctions, that we can generate using PortCullis:

```bash
runPortcullis.sh merged_B73_AB10.bam
```

### 3. Consolidate transcripts and prepare files

Here, we will use all the transcript assemblies generated in the previous step, along with splice junctions and pool them together. The redundant transcripts are also removed. Using consolidated transcripts, we will find full length trasncritps by BLAST searching against SwissProt proteins, and detect ORFs using TransDecoder.

This is all accomplished using this [script](copy-files-for-mikado.sh):

```bash
copy-files-for-mikado.sh B73_Ab10
```

### 4. Finalize evidence-based predictions

As a last step, we will run Mikado to finalize the predictions. It is done using this [script](finalize-files-and-pick.sh):


```bash
finalize-files-and-pick.sh B73_Ab10
```

### 5. Post processing of evidence-based predictions

Kapeel's section
