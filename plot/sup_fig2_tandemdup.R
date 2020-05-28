setwd("/Users/jianingliu/Downloads/genome/")

tandem.repeat <-read.table('v2.repeat.bed')
ngap.final <-read.table('B73Ab10.Ngap.bed')
ngap.pacbio <-read.table('ab10.pacbio.final.ngap.bed')
ngap.ont.contig <-read.table('ont.contig.gap.100kb.bed')
ngap.ont <-read.table('ab10.pacbio.final.ngap.bed')
ont.read.gap <-read.table('B73Ab10.ont.readgap3.bed')
pacbio.read.gap <-read.table('B73Ab10.pacbio.readgap2.bed')
illumina.read.gap <-read.table('B73Ab10.illumina.readgap2.bed')
b73pacbio.read.gap <-read.table('PacBio_error_corrected.3reads.bed')
b73illumina.read.gap <-read.table('b73_illumina.read.gap.bed')
pacbio.readcov <-read.table('PacBio_error_corrected.10kb.bedgraph')
pacbio.readcov[pacbio.readcov$V4>25,]$V4 <- 25
ont.readcov <-read.table('ONT_error_corrected.10kb.bedgraph')
ont.readcov[ont.readcov$V4>25,]$V4 <- 25

tiff("chr8-tandemduplication.tif", width = 8, height = 4, units = 'in', res = 500)
chr3=GRanges(seqnames="chr8",ranges=IRanges(start = 0, end = 182821033))
zoom.region <- toGRanges(data.frame("chr8", 31000000,33500000))
region<-plotKaryotype(chr3,cytobands = Ab10.cytobands, plot.params = pp,zoom=zoom.region)

kpRect(region, data=toGRanges(tandem.repeat),x0=tandem.repeat$V2,x1=tandem.repeat$V3,y0=0, y1=1,lwd=0.07,r0=0.85, r1=0.92,col='#070206',border='#070206')
kpRect(region, data=toGRanges(ngap.final),x0=ngap.final$V2,x1=ngap.final$V3,y0=0, y1=1,lwd=0.25,r0=0.75, r1=0.82,col='#070206',border='#070206')
kpRect(region, data=toGRanges(ngap.pacbio),x0=ngap.pacbio$V2,x1=ngap.pacbio$V3,y0=0, y1=1,lwd=0.25,r0=0.65, r1=0.72,col='#434554',border='#434554')
kpHeatmap(region, data=toGRanges(pacbio.readcov),x0=pacbio.readcov$V2,x1=pacbio.readcov$V3,y=pacbio.readcov$V4,r0=0.55, r1=0.63,colors = c('#FFFFFF','#d4d5dc','#8d8fa3','#71748d','#434554'))
kpRect(region, data=toGRanges(pacbio.read.gap),x0=pacbio.read.gap$V2,x1=pacbio.read.gap$V3,y0=0, y1=1,lwd=0.07,r0=0.45, r1=0.52,col='#434554',border='#434554')
kpRect(region, data=toGRanges(ngap.ont.contig),x0=ngap.ont.contig$V2,x1=ngap.ont.contig$V3,y0=0, y1=1,lwd=0.35,r0=0.35, r1=0.42,col='#620808',border='#620808')
kpHeatmap(region, data=toGRanges(ont.readcov),x0=ont.readcov$V2,x1=ont.readcov$V3,y=ont.readcov$V4,r0=0.25, r1=0.33,colors = c('#FFFFFF','#c09c9c','#915252','#620808','#3a0404'))
#kpRect(region, data=toGRanges(ont.read.gap),x0=ont.read.gap$V2,x1=ont.read.gap$V3,y0=0, y1=1,lwd=0.07,r0=0.45, r1=0.52,col='#667180',border='#667180')
kpRect(region, data=toGRanges(ont.read.gap),x0=ont.read.gap$V2,x1=ont.read.gap$V3,y0=0, y1=1,lwd=0.07,r0=0.15, r1=0.22,col='#620808',border='#620808')
kpRect(region, data=toGRanges(illumina.read.gap),x0=illumina.read.gap$V2,x1=illumina.read.gap$V3,y0=0, y1=1,lwd=0.07,r0=0.05, r1=0.12,col='#003e21',border='#003e21')
#kpRect(region, data=toGRanges(b73pacbio.read.gap),x0=b73pacbio.read.gap$V2,x1=b73pacbio.read.gap$V3,y0=0, y1=1,lwd=0.07,r0=0.15, r1=0.22,col='#8b482d',border='#8b482d')
#kpRect(region, data=toGRanges(b73illumina.read.gap),x0=b73illumina.read.gap$V2,x1=b73illumina.read.gap$V3,y0=0, y1=1,lwd=0.17,r0=0.05, r1=0.12,col='#8b482d',border='#8b482d')
dev.off()

kpAddBaseNumbers(region, tick.dist = 500000, add.units = FALSE,digits=1,cex=1.5,tick.len = 10,minor.tick.dist = 100000, minor.tick.len = 5,clipping=TRUE,cex.axis=3,cex.lab=3,cex.main=3,cex.sub=3)

