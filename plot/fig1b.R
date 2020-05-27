library(karyoploteR)
library(data.table)
library(ggplot2)

setwd("/Users/jianingliu/Downloads/karyoploter/data")
#map9: 57,737,000-61,037,000
#map9:57,443,032-60,752,269
pp <- getDefaultPlotParams(plot.type = 2)
pp$data2height <- 90
pp$ideogramheight <- 0.5
pp$data1inmargin<-0
chr9_ontcontig <- read.table('ont_contig_alignment.bed')[,c('V1','V2','V3')]
chr9_ontcontig$V1 = "chr9"
chr9_pacbiocontig <- read.table('pacbio_contig_alignment.bed')[,c('V1','V2','V3')]
head(chr9_pacbiocontig)
extra_contig<-data.frame(9,57920000,58650000)
names(extra_contig)<-c("V1","V2",'V3')
chr9_pacbiocontig <- rbind(chr9_pacbiocontig, extra_contig)
chr9_pacbiocontig$V1 = "chr9"
chr9_merged.CentC <- read.table('/Users/jianingliu/Downloads/Merged_Final/CentC_1kb.heatmap')
chr9_merged.CentC[chr9_merged.CentC$V4>1,]$V4 <- 1
max(chr9_merged.CentC$V4)
head(chr9_merged.CentC)

tiff("fig1b_a.tiff", width = 8, height = 8, units = 'in', res = 500)
chr9=GRanges(seqnames="chr9",ranges=IRanges(start = 0, end = 161994764))
zoom.region <- toGRanges(data.frame("chr9", 57443032,60752269))
chr9_region<-plotKaryotype(chr9, plot.params = pp,zoom=zoom.region)
kpAddBaseNumbers(chr9_region, tick.dist = 500000, add.units = FALSE,minor.tick.len=0)
kpHeatmap(chr9_region, toGRanges(chr9_merged.CentC), y=chr9_merged.CentC$V4, colors = c(transparent('#FFFFFF', amount=1),'#404246'), r0=0.0, r1=0.15)
kpPlotRegions(chr9_region, data=chr9_pacbiocontig, col="#FFCCAA", layer.margin = 0.05, border="#FFCCAA", r0=0.3, r1=0.6)
dev.off()

tiff("/Users/jianingliu/Downloads/karyoploter/data/fig2b_scale.tiff", width = 8, height = 2, units = 'in', res = 500)
chr9=GRanges(seqnames="chr9",ranges=IRanges(start = 0, end = 161994764))
zoom.region <- toGRanges(data.frame("chr9", 57443032,60752269))
chr9_region<-plotKaryotype(chr9, plot.params = pp,zoom=zoom.region)
kpAddBaseNumbers(chr9_region, tick.dist = 500000, add.units = FALSE,digits=1,cex=1,tick.len = 25,minor.tick.dist = 100000, minor.tick.len = 12)
kpPlotRegions(chr9_region, data=chr9_pacbio.gap,r0=0.15, r1=0.3,col='#0000CD')
kpHeatmap(chr9_region, toGRanges(chr9_pacbio.CentC), y=chr9_pacbio.CentC$V4, colors = c(transparent('#FFFFFF', amount=1),'#404246'), r0=0.0, r1=0.15)
dev.off()


