dat <- data.frame(read.table('repeat_assembly.stats'))
head(dat)
ggplot(dat, aes(x=V1, y=V2/1000000,color=V1)) + geom_boxplot() +scale_color_manual(values=c("#620808", "#8d8fa3", "#211e1a")) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1,size=17),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) +
  xlab("Repeat assembly") + ylab("Repeat array size (Mb)")

