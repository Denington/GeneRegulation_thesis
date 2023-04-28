#### differential RAD21 binding with csaw
#looks much better when normalising for sample specific trended biasss
library(csaw)
library(rtracklayer)

setwd('~/Documents/Projects/Thymocyte_HiC/chip_seqs/bams/')

blacklist <- import.bed('../../Thesis/external_data/mm9_blacklist.bed')
bad.bins <- import.bed('../../bad_bins/CD69negDPWT.bad_HiC_bins.10kb.bed')
blacklist <- c(blacklist,bad.bins)

standard.chr <- paste0("chr", c(1:19))
param <- readParam(minq=25, discard=blacklist,restrict=standard.chr)

bam.files <- c('5_Rad21CD69negDPWTR1_mapped_sorted_RemoveDuplicates.bam','Rad21CD69negDPWTR2_mapped_sorted_RemoveDuplicates.bam',
              'Rad21CD69negDPCTCFKOR1_mapped_sorted_RemoveDuplicates.bam','Rad21CD69negDPCTCFKOR2_mapped_sorted_RemoveDuplicates.bam',
              '6_Rad21CD69pos4SPWTR1_mapped_sorted_RemoveDuplicates.bam','6_Rad21CD69pos4SPWTR2_mapped_sorted_RemoveDuplicates.bam')

Status= c('WT','WT','CTCFKO','CTCFKO','WT4','WT4')
data.frame(BAM=bam.files,  Status= Status)

#####
x <- correlateReads(bam.files, param=reform(param, dedup=TRUE))
frag.len <- which.max(x) - 1
frag.len
#
plot(1:length(x)-1, x, xlab="Delay (bp)", ylab="CCF", type="l")
abline(v=frag.len, col="red")
text(x=frag.len, y=min(x), paste(frag.len, "bp"), pos=4, col="red")
######

win.data <- windowCounts(bam.files, param=param, width=150, ext=frag.len)  #feels a bit short almost- like should be ~500bp?
#150 bp seems like a good distance from trying a couple
#nucleosome ~ 150bp
#and the actually binding profile looks like ~300bp so hmm. But feels like should be a short ~20bp thing
win.data

# 
bins <- windowCounts(bam.files , bin=TRUE, width=2000, param=param) 
filter.stat <- filterWindowsGlobal(win.data, bins)
min.fc <- 3
keep <- filter.stat$filter > log2(min.fc)
summary(keep)
####

filtered.data <- win.data[keep,]

#
win.ab <- scaledAverage(filtered.data)
adjc <- calculateCPM(filtered.data, use.offsets=FALSE)
logfc <- adjc[,4] - adjc[,1]
smoothScatter(win.ab, logfc, ylim=c(-6, 6), xlim=c(0, 5),
              xlab="Average abundance", ylab="Log-fold change")

lfit <- smooth.spline(logfc~win.ab, df=5)
o <- order(win.ab)
lines(win.ab[o], fitted(lfit)[o], col="red", lty=2)

##

filtered.data <- normOffsets(filtered.data)
offsets <- assay(filtered.data, "offset")
head(offsets)

norm.adjc <- calculateCPM(filtered.data, use.offsets=TRUE)
norm.fc <- norm.adjc[,4]-norm.adjc[,1]
smoothScatter(win.ab, norm.fc, ylim=c(-6, 6), xlim=c(0, 5),
              xlab="Average abundance", ylab="Log-fold change")

lfit <- smooth.spline(norm.fc~win.ab, df=5)
lines(win.ab[o], fitted(lfit)[o], col="red", lty=2)
###

Status <- factor(Status)
design <- model.matrix(~0+Status)
colnames(design) <- levels(Status)
design

##

library(edgeR)

y <- asDGEList(filtered.data)
summary(y)
y <- estimateDisp(y, design)
summary(y$trended.dispersion)

plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$df.prior)

plotQLDisp(fit)

#plotMDS(cpm(y, log=TRUE), top=10000, labels=Status,
#        col=c("red", "blue")[as.integer(Status)])
plotMDS(cpm(y, log=TRUE), top=10000, labels=c('CD69nDPWT', 'CD69nDPWT','CD69nDPCTCFKO','CD69nDPCTCFKO','CD69pSPWT','CD69pSPWT'),
        col=c("red", "blue")[as.integer(Status)])
#
contrast <- makeContrasts(CTCFKO-WT, levels=design)
res <- glmQLFTest(fit, contrast=contrast)

plotMD(res,xlab='Log2 CPM',ylab='Log2 Fold Change',hl.cex=0.25,bg.cex=0.1,ylim=c(-8,8),legend=F)
abline(h=0,col='red',lwd=1, lty=2)


merged <- mergeResults(filtered.data, res$table, tol=100, 
                       merge.args=list(max.width=5000))
merged$regions
#
tabcom <- merged$combined
is.sig <- tabcom$FDR <= 0.05
summary(is.sig)

table(tabcom$direction[is.sig])

tabbest <- merged$best
is.sig.pos <- (tabbest$rep.logFC > 0)[is.sig]
summary(is.sig.pos)



#
out.ranges <- merged$regions
mcols(out.ranges) <- DataFrame(tabcom,
                               best.pos=mid(ranges(rowRanges(filtered.data[tabbest$rep.test]))),
                               best.logFC=tabbest$rep.logFC)

out.df <- as.data.frame(out.ranges)
write.table(out.df,'~/Documents/Projects/Thymocyte_HiC/Differential_chip/WTvsCTCFko.RAD21.bed',sep='\t',row.names = F,quote = F)
####s

z$status <- 'Stable'
z$status[out.ranges$FDR < 0.05,] <- 'Dev'
z$status[ctcfkoout.ranges$FDR < 0.05] <- 'CTCFko'
z$status[ctcfkoout.ranges$FDR < 0.05 & out.ranges$FDR < 0.05] <- 'Both'
ggplot(z,aes(x,y,color=status)) + geom_point(alpha=0.3) + ggsci::scale_color_jco()
cor.test(z$x,z$y)

ggplot() + geom_point(data = z[z$status== 'Stable',], aes(x=CTCFko,y=Dev),size=0.25,alpha=0.5) + 
  geom_point(data=z[z$status != 'Stable',] , aes(x=CTCFko,y=dev,color=status),size=0.5,alpha=0.7) + theme_classic() + 
  scale_colour_manual(values = c('red3','blue3',"dodgerblue")) + geom_hline(yintercept = 0,linetype='dashed',size=0.25) + 
  geom_vline(xintercept = 0,linetype='dashed',size=0.25) + coord_cartesian(xlim = c(-6,6),ylim = c(-6,6)) 
  
