#Script for comparing phyloP scores as number of accelerated regions at a cne increases, including and excluding the accelerated species
library(ggplot2)
library(genomation)
library(rtracklayer)

#mypal= ggsci::pal_gsea("default", alpha = 0.7)(9)
setwd('~/Documents/Projects/Accelerated_Regions/phyloP.bws/')

cnes <- import.bed('../CNEs.July17.bed')
cnes_ages <- import.bed('../generalstuff/cne_conservation/CNEs.conserved_to.bed')
cnes_ages <- cnes_ages[unique(queryHits(findOverlaps(cnes_ages,cnes)))]

#######################
#### To start with, metaplots of 100way conservation

#pp100.way.0.w <- ScoreMatrixBin('0accel',
#              resize(cnes_ages,2*width(cnes_age),fix='center'),bin.num = 51,type="bw",strand.aware = F)
#pp100.way.0.w@.Data <- pp100.way.0.w@.Data[rowSums(pp100.way.0.w@.Data) > 1,]

#pp100.way.0.wo <- ScoreMatrixBin('0accelwo',
#                              resize(cnes_ages,2*width(cnes_age),fix='center'),bin.num = 51,type="bw",strand.aware = F)
#pp100.way.0.wo@.Data <- pp100.way.0.wo@.Data[rowSums(pp100.way.0.wo@.Data) > 1,]
#
pp100.way.1.w <- ScoreMatrixBin('1.accels.including.accs.mammal.bw',
                                resize(cnes_ages,2*width(cnes_ages),fix='center'),bin.num = 51,type="bw",strand.aware = F)
pp100.way.1.w@.Data <- pp100.way.1.w@.Data[rowSums(pp100.way.1.w@.Data) > 1,]

pp100.way.1.wo <- ScoreMatrixBin('1.accels.all_but.mammal.bw',
                                 resize(cnes_ages,2*width(cnes_ages),fix='center'),bin.num = 51,type="bw",strand.aware = F)
pp100.way.1.wo@.Data <- pp100.way.1.wo@.Data[rowSums(pp100.way.1.wo@.Data) > 1,]
#
pp100.way.2.w <- ScoreMatrixBin('2.accels.including.accs.mammal.bw',
                                resize(cnes_ages,2*width(cnes_ages),fix='center'),bin.num = 51,type="bw",strand.aware = F)
pp100.way.2.w@.Data <- pp100.way.2.w@.Data[rowSums(pp100.way.2.w@.Data) > 1,]

pp100.way.2.wo <- ScoreMatrixBin('2.accels.all_but.accs.mammal.wig',
                                 resize(cnes_ages,2*width(cnes_ages),fix='center'),bin.num = 51,type="bw",strand.aware = F)
pp100.way.2.wo@.Data <- pp100.way.2.wo@.Data[rowSums(pp100.way.2.wo@.Data) > 1,]
#
pp100.way.3.w <- ScoreMatrixBin('3.accels.including.accs.mammal.bw',
                                resize(cnes_ages,2*width(cnes_ages),fix='center'),bin.num = 51,type="bw",strand.aware = F)
pp100.way.3.w@.Data <- pp100.way.3.w@.Data[rowSums(pp100.way.3.w@.Data) > 1,]

pp100.way.3.wo <- ScoreMatrixBin('3.accels.all_but.mammal.bw',
                                 resize(cnes_ages,2*width(cnes_ages),fix='center'),bin.num = 51,type="bw",strand.aware = F)
pp100.way.3.wo@.Data <- pp100.way.3.wo@.Data[rowSums(pp100.way.3.wo@.Data) > 1,]

#
pp100.way.4.w <- ScoreMatrixBin('4.accels.including.accs.mammal.bw',
                                resize(cnes_ages,2*width(cnes_ages),fix='center'),bin.num = 51,type="bw",strand.aware = F)
pp100.way.4.w@.Data <- pp100.way.4.w@.Data[rowSums(pp100.way.4.w@.Data) > 1,]

pp100.way.4.wo <- ScoreMatrixBin('4.accels.all_but.mammal.bw',
                                 resize(cnes_ages,2*width(cnes_ages),fix='center'),bin.num = 51,type="bw",strand.aware = F)
pp100.way.4.wo@.Data <- pp100.way.4.wo@.Data[rowSums(pp100.way.4.wo@.Data) > 1,]

#
pp100.way.5.w <- ScoreMatrixBin('5.accels.including.accs.mammal.bw',
                                resize(cnes_ages,2*width(cnes_ages),fix='center'),bin.num = 51,type="bw",strand.aware = F)
pp100.way.5.w@.Data <- pp100.way.5.w@.Data[rowSums(pp100.way.5.w@.Data) > 1,]

pp100.way.5.wo <- ScoreMatrixBin('5.accels.all_but.mammal.bw',
                                 resize(cnes_ages,2*width(cnes_ages),fix='center'),bin.num = 51,type="bw",strand.aware = F)
pp100.way.5.wo@.Data <- pp100.way.5.wo@.Data[rowSums(pp100.way.5.wo@.Data) > 1,]


#

pp100.way.6.w <- ScoreMatrixBin('6.accels.including.accs.mammal.bw',
                                resize(cnes_ages,2*width(cnes_ages),fix='center'),bin.num = 51,type="bw",strand.aware = F)
pp100.way.6.w@.Data <- pp100.way.6.w@.Data[rowSums(pp100.way.6.w@.Data) > 1,]

pp100.way.6.wo <- ScoreMatrixBin('6.accels.all_but.mammal.bw',
                                 resize(cnes_ages,2*width(cnes_ages),fix='center'),bin.num = 51,type="bw",strand.aware = F)
pp100.way.6.wo@.Data <- pp100.way.6.wo@.Data[rowSums(pp100.way.6.wo@.Data) > 1,]


####

pp100.way.7.w <- ScoreMatrixBin('7.accels.including.accs.mammal.bw',
                                resize(cnes_ages,2*width(cnes_ages),fix='center'),bin.num = 51,type="bw",strand.aware = F)
pp100.way.7.w@.Data <- pp100.way.7.w@.Data[rowSums(pp100.way.7.w@.Data) > 1,]

pp100.way.7.wo <- ScoreMatrixBin('7.accels.all_but.mammal.bw',
                                 resize(cnes_ages,2*width(cnes_ages),fix='center'),bin.num = 51,type="bw",strand.aware = F)
pp100.way.7.wo@.Data <- pp100.way.7.wo@.Data[rowSums(pp100.way.7.wo@.Data) > 1,]


##########################
oneAcc <- data.frame(scores = c(rowMeans(pp100.way.1.w@.Data),rowMeans(pp100.way.1.wo@.Data)), Accels=c(rep('1',dim(pp100.way.1.w)[1]), rep('1',dim(pp100.way.1.wo)[1])),
           AccInc=c(rep('With',dim(pp100.way.1.w)[1]), rep('Without',dim(pp100.way.1.wo)[1])))

twoAcc <-data.frame(scores = c(rowMeans(pp100.way.2.w@.Data),rowMeans(pp100.way.2.wo@.Data)), Accels=c(rep('2',dim(pp100.way.2.w)[1]), rep('2',dim(pp100.way.2.wo)[1])),
                     AccInc=c(rep('With',dim(pp100.way.2.w)[1]), rep('Without',dim(pp100.way.2.wo)[1])))

threeAcc <-data.frame(scores = c(rowMeans(pp100.way.3.w@.Data),rowMeans(pp100.way.3.wo@.Data)), Accels=c(rep('3',dim(pp100.way.3.w)[1]), rep('3',dim(pp100.way.3.wo)[1])),
                     AccInc=c(rep('With',dim(pp100.way.3.w)[1]), rep('Without',dim(pp100.way.3.wo)[1])))

fourAcc <-data.frame(scores = c(rowMeans(pp100.way.4.w@.Data),rowMeans(pp100.way.4.wo@.Data)), Accels=c(rep('4',dim(pp100.way.4.w)[1]), rep('4',dim(pp100.way.4.wo)[1])),
                    AccInc=c(rep('With',dim(pp100.way.4.w)[1]), rep('Without',dim(pp100.way.4.wo)[1])))

fiveAcc <-data.frame(scores = c(rowMeans(pp100.way.5.w@.Data),rowMeans(pp100.way.5.wo@.Data)), Accels=c(rep('5',dim(pp100.way.5.w)[1]), rep('5',dim(pp100.way.5.wo)[1])),
                    AccInc=c(rep('With',dim(pp100.way.5.w)[1]), rep('Without',dim(pp100.way.5.wo)[1])))


sixAcc <-data.frame(scores = c(rowMeans(pp100.way.6.w@.Data),rowMeans(pp100.way.6.wo@.Data)), Accels=c(rep('6',dim(pp100.way.6.w)[1]), rep('6',dim(pp100.way.6.wo)[1])),
                    AccInc=c(rep('With',dim(pp100.way.6.w)[1]), rep('Without',dim(pp100.way.6.wo)[1])))
sevenAcc <- data.frame(scores = c(rowMeans(pp100.way.7.w@.Data),rowMeans(pp100.way.7.wo@.Data)), Accels=c(rep('7',dim(pp100.way.7.w)[1]), rep('7',dim(pp100.way.7.wo)[1])),
                         AccInc=c(rep('With',dim(pp100.way.7.w)[1]), rep('Without',dim(pp100.way.7.wo)[1])))

ggplot(rbind(oneAcc,twoAcc,threeAcc,fourAcc,fiveAcc,sixAcc,sevenAcc),aes(Accels,scores,fill=AccInc)) + geom_boxplot() + ggpubr::stat_compare_means()
###########################

testgen <- ScoreMatrixBin('~/Documents/Projects/Accelerated_Regions/1.accels.excluding.woofboy.bw',
                                  resize(cnes_ages[cnes_ages$score == 3],2*width(cnes_ages[cnes_ages$score == 3]),fix='center'),bin.num = 51,type="bw",strand.aware = F)

x <-  cnes_ages[cnes_ages$score == 3]

testgen <- ScoreMatrixBin('~/Documents/Projects/Accelerated_Regions/1.accels.excluding.woofboy.bw',
                         resize(x[order(width(x),decreasing = T)],1500,fix='center'),bin.num = 51,type="bw",strand.aware = F)
heatMatrix(testgen,winsorize = c(1,99)) #,col=mypal)