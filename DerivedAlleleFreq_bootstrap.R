################################################
#File Name: DerivedAlleleFreq_bootstrap.R
#Usage: Do bootstrap of the DAF and Draw plot
#Example of length DAF
################################################

#!/bin/env Rscript

#input data with constant format
data <- read.table('DAF_byLen.bed3+',
                   col.names = c("chrom", "start","end", "SVname", 
                                 "ancestorAlleleCount",  "derivedAlleleCount", "derivedAlleleFreq", "type"))
data$TotalN <- data$ancestorAlleleCount + data$derivedAlleleCount

#determine the sampling count for bootstrap to contain most samples
quantile(data$TotalN, probs = c(0.001, 0.05, 0.1, 0.2))

#count numbers in each bin region
sfs <- function(df, br = seq(0, 1, 0.05) * C) {
  x <- hist(df, breaks = br, plot = FALSE)
  return(x$counts)
}
#shuffle with replace
random<-function(vctor,C){
  rst<-c()
  for (i in (1:length(vctor))){
    bsample<-sample(c(rep(1,vctor[i]),rep(0,C-vctor[i])),replace=T)
    rst<-c(rst,sum(bsample))
  }
  return(rst)
}
lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}

C=500
data1 <- subset(data, TotalN >= C)
tmp <- as.matrix(cbind(data1$ancestorAlleleCount, data1$derivedAlleleCount))

#rescue situations if derived allele count=1
subset <- apply(tmp, 1, function(t){rhyper(10, t[2]-2, t[1], C-2)})
subset = subset + 2
subset[is.na(subset)]<-1
i = 1
split <- split(subset[i, ], list(data1$type))
below1kb<-list()
over1kb<-list()

#why sample！=C：because alt count!=0，when ref！=anc，ref=derived！=C
for(i in 1:1000){
  sample<-random(split$below1kb,C=C)
  below1kb<-lappend(below1kb,sample[which(sample!=0 & sample!=C)])
  sample<-random(split$over1kb,C=C)
  over1kb<-lappend(over1kb,sample[which(sample!=0 & sample!=C)])
}

SFS_below1kb<- sapply(below1kb,sfs)
SFSd_below1kb<-apply(SFS_below1kb, 2, function(n) n/sum(n))
SFS_over1kb<-sapply(over1kb,sfs)
SFSd_over1kb<- apply(SFS_over1kb, 2, function(n) n/sum(n))

library(ggplot2)
bin<-c('0.0-0.05','0.05-0.1','0.1-0.15','0.15-0.2','0.2-0.25','0.25-0.3','0.3-0.35','0.35-0.4',
       '0.4-0.45','0.45-0.5','0.5-0.55','0.55-0.6','0.6-0.65','0.65-0.7','0.7-0.75','0.75-0.8',
       '0.8-0.85','0.85-0.9','0.9-0.95','0.95-1')

data <- as.data.frame(c(apply(SFSd_below1kb,1,mean),apply(SFSd_over1kb,1,mean)))
colnames(data)<-c('freq')
data$bin<-c(bin,bin)
data$Class<-c(rep('below1kb',length(bin)),rep('over1kb',length(bin)))
data$sd<-c(apply(SFSd_below1kb,1,sd),apply(SFSd_over1kb,1,sd))
pdf("inversion_by_length_DAF.pdf", width = 9, height = 6)
data$Class<-factor(data$Class, levels = c("over1kb","below1kb"))
ggplot(data,aes(x=bin,y=freq,fill=Class)) +
  geom_bar(stat="identity",position="dodge") +
  geom_errorbar(aes(ymax=freq+sd,ymin=freq-sd),position=position_dodge(0.9),width=0.2,size=0.4) +
  scale_fill_brewer(palette="Set1") + theme_bw() +
  theme(axis.text.x = element_text(angle=20, hjust=1))+
  scale_y_continuous(breaks = seq(0,1,0.1),limits=c(0,0.75))+
  theme(legend.position=c(0.8,0.67),text=element_text(size=15),panel.grid =element_blank()) +
  ggtitle("By length")+xlab('Derived Allele Frequency') +
  ylab('Fraction of inversions')
dev.off()
