n50 <- function(len) {
  len.sorted <- rev(sort(len))
  N50 <- len.sorted[cumsum(as.numeric(len.sorted)) >= sum(len.sorted)*0.5][1]
  return(N50)
}

dat <- data.frame(read.table('chr8.10kb.readlength'))

head(dat)
options(scipen=1000000)

ont<-dat[dat$V1=='ONT',]$V2
pacbio<-dat[dat$V1=='PacBio',]$V2
n50(ont)
n50(pacbio)

wilcox.test(ont,pacbio) # where y and x are numeric

tiff("chr8.tif", width = 4, height = 4, units = 'in', res = 500)
ggplot(dat, aes(x=V2/1000, color=V1)) +
  geom_histogram(fill="white", alpha=0.1, position="identity",bins=100) +
  scale_color_manual(values=c( "#620808", "#8d8fa3")) +
  theme_minimal() + geom_density(aes(y=..count../100)) +scale_x_continuous(trans='log10') 
  + scale_y_log10() +
  theme(axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14))  +
  geom_vline(xintercept = n50(ont)/1000,colour = "#999999", linetype="dashed")+
  geom_vline(xintercept = n50(pacbio)/1000,colour = "#E69F00", linetype="dashed") +scale_x_continuous(trans='log10') +
  geom_histogram(fill="white", alpha=0.3, position="identity",bins=100,aes(y=..density..)) +   geom_density(linetype="dashed")
dev.off()

