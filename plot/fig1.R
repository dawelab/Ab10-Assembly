library(karyoploteR)
library(data.table)
library(ggplot2)

setwd("~/karyoploter/data")

pp <- getDefaultPlotParams(plot.type = 2)
pp$data2height <- 90
pp$ideogramheight <- 0.5
pp$data1inmargin<-0
Ab10.cytobands <- toGRanges("Ab10_cytobands.txt")
Ab10.genome <- toGRanges(read.table('chrs.sizes'))
Ab10.contigs <- read.table('chrs_contigs.bed')
full_chrs <- Ab10.contigs[Ab10.contigs$V1 =='chr3'|Ab10.contigs$V1=='chr9',]
Ab10.contigs<-Ab10.contigs[Ab10.contigs$V1!='chr3'&Ab10.contigs$V1!='chr9',]
slice1 <- toGRanges(Ab10.contigs[seq(1, nrow(Ab10.contigs), by = 3),])
slice2 <- toGRanges(Ab10.contigs[seq(2, nrow(Ab10.contigs), by = 3),])
slice3 <- toGRanges(Ab10.contigs[seq(3, nrow(Ab10.contigs), by = 3),])

# read heatmap files
Ab10.Copia.heatmap <- read.table('LTR_Copia_10kb.heatmap')
Ab10.Copia.heatmap.range <- toGRanges(Ab10.Copia.heatmap[,c('V1','V2','V3')])
Ab10.Gypsy.heatmap <- read.table('LTR_Gypsy_10kb.heatmap')
Ab10.Gypsy.heatmap.range <- toGRanges(Ab10.Gypsy.heatmap[,c('V1','V2','V3')])
Ab10.cinful_zeon.heatmap <- read.table('cinful_zeon_10kb.heatmap')
Ab10.cinful_zeon.heatmap.range <- toGRanges(Ab10.cinful_zeon.heatmap[,c('V1','V2','V3')])
Ab10.gene <- read.table('B73_AB10_evd_10kb.heatmap',header = FALSE)
Ab10.gene.range <-toGRanges(Ab10.gene[,c('V1','V2','V3')])

# read repeat and TE bed files 
Ab10.TR1 <- read.table('Filtered_TR-1_total.bed')
Ab10.kno180 <- read.table('Filtered_knob180_total.bed')
Ab10.CentC <- read.table('Filtered_CentC_total.bed')

Ab10.huck <- read.table('huck.bed')
Ab10.opieji <- read.table('opie-ji-super.bed')
Ab10.cinful_zeon <- read.table('cinful_zeon.bed')
Ab10.prem1 <- read.table('prem1-super.bed')
Ab10.grande <- read.table('grande.bed')
Ab10.crm <- read.table('crm.bed')
Ab10.active.centromere <- read.csv('active_centromere.csv',header = TRUE)

#sort(Ab10.gene$V4,decreasing = TRUE)
#quantile(Ab10.gene$V4,probs = c(0.05, 0.85))
Ab10.gene$V4[Ab10.gene$V4 >= 1000] <- 1000
Ab10.Gypsy.heatmap.range <- toGRanges(Ab10.Gypsy.heatmap[,c('V1','V2','V3')])
#sort(Ab10.Gypsy.heatmap$V4,decreasing = TRUE)
#quantile(Ab10.Gypsy.heatmap$V4,probs = c(0.05, 0.85))
Ab10.Gypsy.heatmap$V4[Ab10.Gypsy.heatmap$V4 >= 4000] <- 4000
Ab10.cinful_zeon.heatmap$V4[Ab10.cinful_zeon.heatmap$V4 >= 2000] <- 2000
#quantile(Ab10.cinful_zeon.heatmap$V4,probs = c(0.05, 0.85))

tiff("Fig1_modified1_cifulzeon.tif", width = 10, height = 8, units = 'in', res = 400)

Ab10.cytobands <- toGRanges("Ab10_cytobands.txt")
genome_plot <- plotKaryotype(genome = Ab10.genome,cytobands = Ab10.cytobands, plot.type=2, plot.params = pp,ideogram.plotter = NULL)
kpRect(genome_plot, chr=Ab10.active.centromere$Chr, x0=Ab10.active.centromere$Start, x1=Ab10.active.centromere$End, y0= 0,y1=1, r0=-0.03, r1=0.15,lwd=0.7,border='#FF8C00',col='#FF8C00')
kpHeatmap(genome_plot, Ab10.Gypsy.heatmap.range, y=Ab10.Gypsy.heatmap$V4, colors = c(transparent('#FFFFFF', amount=0),'#556272','#3b444f','#22272d'), r0=0.25, r1=0.5)
kpHeatmap(genome_plot, Ab10.cinful_zeon.heatmap.range, y=Ab10.cinful_zeon.heatmap$V4, colors = c(transparent('#FFFFFF', amount=0),'#401b89','#2e1362','#1b0b3b'), r0=0.5, r1=0.75)
kpHeatmap(genome_plot, Ab10.gene.range, y=Ab10.gene$V4, colors = c(transparent('#FFFFFF', amount=0),'#206e78','#19565d','#154950','#123d43'), r0=0.75, r1=1.0)

kpPlotRegions(genome_plot, data=slice2, col="#ECCBAE", border="#363849",r0=-0.2, r1=-0.1,avoid.overlapping=FALSE)
kpPlotRegions(genome_plot, data=slice1, col="#ECCBAE", border="#363849",r0=-0.3, r1=-0.2,avoid.overlapping=FALSE)
kpPlotRegions(genome_plot, data=slice3, col="#ECCBAE", border="#363849",r0=-0.4, r1=-0.3,avoid.overlapping=FALSE)
kpPlotRegions(genome_plot, data=full_chrs, col="#046C9A", border="#363849",r0=-0.2, r1=-0.1,avoid.overlapping=FALSE)

dev.off()

# plot repeat array (inset)

Ab10.gene <- read.table('B73_AB10_evd.bed',header = FALSE)
tiff("ab10_chr4_TEfamilies.tif", width = 16, height = 8, units = 'in', res = 500)
chr1=GRanges(seqnames="chr4",ranges=IRanges(start = 0, end = 249748697))
zoom.region <- toGRanges(data.frame("chr4", 230000000,235000000))

chr1_region<-plotKaryotype(chr1,cytobands = Ab10.cytobands, plot.params = pp,zoom=zoom.region)

kpRect(chr1_region, data=toGRanges(Ab10.kno180),x0=Ab10.kno180$V2,x1=Ab10.kno180$V3,y0=0, y1=1,lwd=0.00,r0=-0.1, r1=0.12,col='#0000CD',border=NA)
kpRect(chr1_region, data=toGRanges(Ab10.TR1),x0=Ab10.TR1$V2,x1=Ab10.TR1$V3,y0=0, y1=1,lwd=0.00,r0=-0.1, r1=0.12,col='#B22222',border=NA)
kpRect(chr1_region, data=toGRanges(Ab10.huck),x0=Ab10.huck$V2,x1=Ab10.huck$V3,y0=0, y1=1,lwd=0.00,r0=0.24, r1=0.36,col='#7c6e81',border=NA)
kpRect(chr1_region, data=toGRanges(Ab10.cinful_zeon),x0=Ab10.cinful_zeon$V2,x1=Ab10.cinful_zeon$V3,y0=0, y1=1,lwd=0.00,r0=0.12, r1=0.24,col='#667180',border=NA)
kpRect(chr1_region, data=toGRanges(Ab10.opieji),x0=Ab10.opieji$V2,x1=Ab10.opieji$V3,y0=0, y1=1,lwd=0.00,r0=0.36, r1=0.48,col='#85865f',border=NA)
kpRect(chr1_region, data=toGRanges(Ab10.grande),x0=Ab10.grande$V2,x1=Ab10.grande$V3,y0=0, y1=1,lwd=0.00,r0=0.48, r1=0.6,col='#be8b2f',border=NA)
kpRect(chr1_region, data=toGRanges(Ab10.prem1),x0=Ab10.prem1$V2,x1=Ab10.prem1$V3,y0=0, y1=1,lwd=0.00,r0=0.6, r1=0.72,col='#8b482d',border=NA)
kpRect(chr1_region, data=toGRanges(Ab10.gene),x0=Ab10.gene$V2,x1=Ab10.gene$V3,y0=0, y1=1,lwd=0.00,r0=0.72, r1=0.84,col='#26705e',border=NA)

kpAddBaseNumbers(chr1_region, tick.dist = 500000, add.units = FALSE,digits=1,cex=1.5,tick.len = 10,minor.tick.dist = 100000, minor.tick.len = 5,clipping=TRUE)

kpAddLabels(chr1_region, labels="cinful-zeon", r0=0.12, r1=0.24,side="left",label.margin=0.005,cex=1.8)
kpAddLabels(chr1_region, labels="huck", r0=0.24, r1=0.36,side="left",label.margin=0.005,cex=1.8)
kpAddLabels(chr1_region, labels="opie-ji", r0=0.36, r1=0.48,side="left",label.margin=0.005,cex=1.8)
kpAddLabels(chr1_region, labels="grande", r0=0.48, r1=0.6,side="left",label.margin=0.005,cex=1.8)
kpAddLabels(chr1_region, labels="prem1", r0=0.6, r1=0.72,side="left",label.margin=0.005,cex=1.8)
kpAddLabels(chr1_region, labels="gene", r0=0.72, r1=0.84,side="left",label.margin=0.005,cex=1.8)

dev.off()

tiff("ab10_chr7_TEfamilies.tif", width = 16, height = 8, units = 'in', res = 500)
chr1=GRanges(seqnames="chr7",ranges=IRanges(start = 0, end = 184241495))
zoom.region <- toGRanges(data.frame("chr7", 154000000,159000000))
chr1_region<-plotKaryotype(chr1,cytobands = Ab10.cytobands, plot.params = pp,zoom=zoom.region)

kpRect(chr1_region, data=toGRanges(Ab10.kno180),x0=Ab10.kno180$V2,x1=Ab10.kno180$V3,y0=0, y1=1,lwd=0.00,r0=-0.1, r1=0.12,col='#0000CD',border=NA)
kpRect(chr1_region, data=toGRanges(Ab10.TR1),x0=Ab10.TR1$V2,x1=Ab10.TR1$V3,y0=0, y1=1,lwd=0.00,r0=-0.1, r1=0.12,col='#B22222',border=NA)
kpRect(chr1_region, data=toGRanges(Ab10.huck),x0=Ab10.huck$V2,x1=Ab10.huck$V3,y0=0, y1=1,lwd=0.00,r0=0.24, r1=0.36,col='#7c6e81',border=NA)
kpRect(chr1_region, data=toGRanges(Ab10.cinful_zeon),x0=Ab10.cinful_zeon$V2,x1=Ab10.cinful_zeon$V3,y0=0, y1=1,lwd=0.00,r0=0.12, r1=0.24,col='#667180',border=NA)
kpRect(chr1_region, data=toGRanges(Ab10.opieji),x0=Ab10.opieji$V2,x1=Ab10.opieji$V3,y0=0, y1=1,lwd=0.00,r0=0.36, r1=0.48,col='#85865f',border=NA)
kpRect(chr1_region, data=toGRanges(Ab10.grande),x0=Ab10.grande$V2,x1=Ab10.grande$V3,y0=0, y1=1,lwd=0.00,r0=0.48, r1=0.6,col='#be8b2f',border=NA)
kpRect(chr1_region, data=toGRanges(Ab10.prem1),x0=Ab10.prem1$V2,x1=Ab10.prem1$V3,y0=0, y1=1,lwd=0.00,r0=0.6, r1=0.72,col='#8b482d',border=NA)
kpRect(chr1_region, data=toGRanges(Ab10.gene),x0=Ab10.gene$V2,x1=Ab10.gene$V3,y0=0, y1=1,lwd=0.00,r0=0.72, r1=0.84,col='#26705e',border=NA)

kpAddBaseNumbers(chr1_region, tick.dist = 500000, add.units = FALSE,digits=1,cex=1.5,tick.len = 10,minor.tick.dist = 100000, minor.tick.len = 5,clipping=TRUE)

kpAddLabels(chr1_region, labels="cinful-zeon", r0=0.12, r1=0.24,side="left",label.margin=0.005,cex=1.8)
kpAddLabels(chr1_region, labels="huck", r0=0.24, r1=0.36,side="left",label.margin=0.005,cex=1.8)
kpAddLabels(chr1_region, labels="opie-ji", r0=0.36, r1=0.48,side="left",label.margin=0.005,cex=1.8)
kpAddLabels(chr1_region, labels="grande", r0=0.48, r1=0.6,side="left",label.margin=0.005,cex=1.8)
kpAddLabels(chr1_region, labels="prem1", r0=0.6, r1=0.72,side="left",label.margin=0.005,cex=1.8)
kpAddLabels(chr1_region, labels="gene", r0=0.72, r1=0.84,side="left",label.margin=0.005,cex=1.8)

dev.off()

