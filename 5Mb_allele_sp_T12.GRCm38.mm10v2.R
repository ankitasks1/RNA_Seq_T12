#RNA-Seq based genotyping
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb")
#bedtools makewindows -g /home/ankitv/ref_av/mm10/mm10.chrom.sizes -w 5000000 > mm10_5Mb.txt
#Calculate count coverage before split
#bedtools bamtobed -i JB1_WT_Rep1_Aligned.sortedByCoord.out.bam > JB1_WT_Rep1_Aligned.sortedByCoord.out.bed 

#bedtools bamtobed -i JB1_WT_Rep2_Aligned.sortedByCoord.out.bam > JB1_WT_Rep2_Aligned.sortedByCoord.out.bed

#bedtools bamtobed -i JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bam > JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bed


#bedtools bamtobed -i JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bam > JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bed

#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' JB1_WT_Rep1_Aligned.sortedByCoord.out.bed > JB1_WT_Rep1_Aligned.sortedByCoord.out_chr.bed

#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' JB1_WT_Rep2_Aligned.sortedByCoord.out.bed > JB1_WT_Rep2_Aligned.sortedByCoord.out_chr.bed

#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bed > JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out_chr.bed

#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bed > JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out_chr.bed


#bedtools coverage -a mm10_5Mb.txt -b JB1_WT_Rep1_Aligned.sortedByCoord.out_chr.bed > JB1_WT_Rep1_Aligned.sortedByCoord.out_chr.cov.bed
#bedtools coverage -a mm10_5Mb.txt -b JB1_WT_Rep2_Aligned.sortedByCoord.out_chr.bed > JB1_WT_Rep2_Aligned.sortedByCoord.out_chr.cov.bed
#bedtools coverage -a mm10_5Mb.txt -b JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out_chr.bed > JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out_chr.cov.bed
#bedtools coverage -a mm10_5Mb.txt -b JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out_chr.bed > JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out_chr.cov.bed

#paste JB1_WT_Rep1_Aligned.sortedByCoord.out_chr.cov.bed JB1_WT_Rep2_Aligned.sortedByCoord.out_chr.cov.bed JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out_chr.cov.bed JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out_chr.cov.bed > mm10_5mb_bsplit.txt

bcovdata5 <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/beforesplit_featurecount/bam/mm10_5mb_bsplit.txt", header=FALSE)
head(bcovdata5)
tail(bcovdata5)
dim(bcovdata5)

bcovdata5_reaaranged <- bcovdata5[,c(1:3, 4, 11, 18, 25)]
head(bcovdata5_reaaranged)
colnames(bcovdata5_reaaranged) <- c("chr","start","end","WT1", "WT2", "ZFP57KO1", "ZFP57KO2")
head(bcovdata5_reaaranged)
colSums(bcovdata5_reaaranged[,4:7])



#     WT1      WT2 ZFP57KO1 ZFP57KO2 
#23729676 24283291 24435673 24831197


#Bulk (Actual Bulk : Bulk_mm10)
#bedtools bamtobed -i JB1_WT_Rep1_Aligned.sortedByCoord.out.bam > JB1_WT_Rep1_Aligned.sortedByCoord.out.bed
#bedtools bamtobed -i JB1_WT_Rep2_Aligned.sortedByCoord.out.bam > JB1_WT_Rep2_Aligned.sortedByCoord.out.bed
#bedtools bamtobed -i JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bam > JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bed
#bedtools bamtobed -i JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bam > JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bed

#bedtools coverage -a mm10_5Mb.txt -b JB1_WT_Rep1_Aligned.sortedByCoord.out.bed > JB1_WT_Rep1.cov_Aligned.sortedByCoord.out.bed
#bedtools coverage -a mm10_5Mb.txt -b JB1_WT_Rep2_Aligned.sortedByCoord.out.bed > JB1_WT_Rep2.cov_Aligned.sortedByCoord.out.bed
#bedtools coverage -a mm10_5Mb.txt -b JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bed > JB1_ZFP57_KO_Rep1.cov_Aligned.sortedByCoord.out.bed
#bedtools coverage -a mm10_5Mb.txt -b JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bed > JB1_ZFP57_KO_Rep2.cov_Aligned.sortedByCoord.out.bed

#paste JB1_WT_Rep1.cov_Aligned.sortedByCoord.out.bed JB1_WT_Rep2.cov_Aligned.sortedByCoord.out.bed JB1_ZFP57_KO_Rep1.cov_Aligned.sortedByCoord.out.bed JB1_ZFP57_KO_Rep2.cov_Aligned.sortedByCoord.out.bed > mm10_5mb_bulksplit.txt

setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2")
bulkcovdata5 <- read.table("mm10_5mb_bulksplit.txt", header=FALSE)
head(bulkcovdata5)
tail(bulkcovdata5)
dim(bulkcovdata5)
bulkcovdata5_reaaranged <- bulkcovdata5[,c(1:3, 4, 11, 18, 25)]
head(bulkcovdata5_reaaranged)
colnames(bulkcovdata5_reaaranged) <- c("chr","start","end","JB1_WT_Rep1","JB1_WT_Rep2","JB1_ZFP57_KO_Rep1","JB1_ZFP57_KO_Rep2")
head(bulkcovdata5_reaaranged)
colSums(bulkcovdata5_reaaranged[,4:7])
#WT1 WT2 ZFP57KO1 ZFP57KO2
#25189744          25429072          25605580          26043468



#bedtools bamtobed -i JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bam > JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bed
#bedtools bamtobed -i JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bam > JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bed
#bedtools bamtobed -i JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam > JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bed
#bedtools bamtobed -i JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam > JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bed
#bedtools bamtobed -i JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bam > JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bed
#bedtools bamtobed -i JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bam > JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bed
#bedtools bamtobed -i JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam > JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bed
#bedtools bamtobed -i JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam > JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bed

#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bed > JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.B6_chr.bed
#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bed > JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1_chr.bed
#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bed > JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.B6_chr.bed
#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bed > JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1_chr.bed
#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bed > JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.B6_chr.bed
#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bed > JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1_chr.bed
#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bed > JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.B6_chr.bed
#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bed > JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1_chr.bed


#bedtools coverage -a mm10_5Mb.txt -b JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.B6_chr.bed > JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.B6_chr.cov.bed
#bedtools coverage -a mm10_5Mb.txt -b JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1_chr.bed > JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1_chr.cov.bed
#bedtools coverage -a mm10_5Mb.txt -b JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.B6_chr.bed > JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.B6_chr.cov.bed
#bedtools coverage -a mm10_5Mb.txt -b JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1_chr.bed > JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1_chr.cov.bed
#bedtools coverage -a mm10_5Mb.txt -b JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.B6_chr.bed > JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.B6_chr.cov.bed
#bedtools coverage -a mm10_5Mb.txt -b JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1_chr.bed > JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1_chr.cov.bed
#bedtools coverage -a mm10_5Mb.txt -b JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.B6_chr.bed > JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.B6_chr.cov.bed
#bedtools coverage -a mm10_5Mb.txt -b JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1_chr.bed > JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1_chr.cov.bed

#paste JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.B6_chr.cov.bed JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.B6_chr.cov.bed JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1_chr.cov.bed JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1_chr.cov.bed JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.B6_chr.cov.bed JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.B6_chr.cov.bed JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1_chr.cov.bed JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1_chr.cov.bed  > mm10_5mb_covdata5.txt


setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb")
covdata5 <- read.table("mm10_5mb_covdata5.txt", header=FALSE)
head(covdata5)
tail(covdata5)
dim(covdata5)
covdata5_reaaranged <- covdata5[,c(1:3, 4, 11, 18, 25, 32, 39, 46, 53)]
head(covdata5_reaaranged)
covdata5_reaaranged["id"] <- data.frame(paste(covdata5_reaaranged$V1, 
                                              covdata5_reaaranged$V2, 
                                              covdata5_reaaranged$V3, sep =  "%"))
head(covdata5_reaaranged)
covdata5_reaaranged <- covdata5_reaaranged[,c(12,4:11)]
rownames(covdata5_reaaranged) <- covdata5_reaaranged[,1]
covdata5_reaaranged <- covdata5_reaaranged[,-1]
colnames(covdata5_reaaranged) <- c("JB1_WT_Rep1.B6","JB1_WT_Rep2.B6","JB1_WT_Rep1.JF1","JB1_WT_Rep2.JF1","JB1_ZFP57_KO_Rep1.B6","JB1_ZFP57_KO_Rep2.B6","JB1_ZFP57_KO_Rep1.JF1","JB1_ZFP57_KO_Rep2.JF1")
head(covdata5_reaaranged)
#colSums(covdata5_reaaranged)
#Factors are the colSums(bulkcovdata5_reaaranged[,4:7]): 25189744          25429072          25605580          26043468

covdata5_reaarangedCPM <- data.frame(cbind(data.frame(rownames(covdata5_reaaranged)),
                                           (covdata5_reaaranged$JB1_WT_Rep1.B6 *1000000)/25189744,
                                           (covdata5_reaaranged$JB1_WT_Rep2.B6*1000000)/25429072,
                                           (covdata5_reaaranged$JB1_WT_Rep1.JF1*1000000)/25189744,
                                           (covdata5_reaaranged$JB1_WT_Rep2.JF1*1000000)/25429072,
                                           (covdata5_reaaranged$JB1_ZFP57_KO_Rep1.B6*1000000)/25605580,
                                           (covdata5_reaaranged$JB1_ZFP57_KO_Rep2.B6*1000000)/26043468,
                                           (covdata5_reaaranged$JB1_ZFP57_KO_Rep1.JF1*1000000)/25605580,
                                           (covdata5_reaaranged$JB1_ZFP57_KO_Rep2.JF1*1000000)/26043468))

head(covdata5_reaarangedCPM)
rownames(covdata5_reaarangedCPM) <- covdata5_reaarangedCPM[,1]
covdata5_reaarangedCPM <- covdata5_reaarangedCPM[,-1]
colnames(covdata5_reaarangedCPM) <- c("JB1_WT_Rep1.B6","JB1_WT_Rep2.B6","JB1_WT_Rep1.JF1","JB1_WT_Rep2.JF1","JB1_ZFP57_KO_Rep1.B6","JB1_ZFP57_KO_Rep2.B6","JB1_ZFP57_KO_Rep1.JF1","JB1_ZFP57_KO_Rep2.JF1")
head(covdata5_reaarangedCPM)
covdata5_reaarangedCPM_sorted <- covdata5_reaarangedCPM[order(rownames(covdata5_reaarangedCPM)),]
head(covdata5_reaarangedCPM_sorted)
covdata5_reaarangedCPM_sortedavg <- data.frame(cbind(data.frame(rownames(covdata5_reaarangedCPM_sorted)),
                                                     (covdata5_reaarangedCPM_sorted$JB1_WT_Rep1.B6 + covdata5_reaarangedCPM_sorted$JB1_WT_Rep2.B6)/2,
                                                     (covdata5_reaarangedCPM_sorted$JB1_WT_Rep1.JF1 + covdata5_reaarangedCPM_sorted$JB1_WT_Rep2.JF1)/2,
                                                     (covdata5_reaarangedCPM_sorted$JB1_ZFP57_KO_Rep1.B6 + covdata5_reaarangedCPM_sorted$JB1_ZFP57_KO_Rep2.B6)/2,
                                                     (covdata5_reaarangedCPM_sorted$JB1_ZFP57_KO_Rep1.JF1 + covdata5_reaarangedCPM_sorted$JB1_ZFP57_KO_Rep2.JF1)/2))

head(covdata5_reaarangedCPM_sortedavg)

rownames(covdata5_reaarangedCPM_sortedavg) <- covdata5_reaarangedCPM_sortedavg[,1]

covdata5_reaarangedCPM_sortedavg <- covdata5_reaarangedCPM_sortedavg[,-1]
colnames(covdata5_reaarangedCPM_sortedavg) <- c("WT_JB1B6","WT_JB1JF1","ZFP57KO_JB1B6","ZFP57KO_JB1JF1")
head(covdata5_reaarangedCPM_sortedavg)
covdata5_reaarangedCPM_sortedavg["matRatioWT"] <- covdata5_reaarangedCPM_sortedavg$WT_JB1JF1/(covdata5_reaarangedCPM_sortedavg$WT_JB1B6 + covdata5_reaarangedCPM_sortedavg$WT_JB1JF1)
covdata5_reaarangedCPM_sortedavg["matRatioZFP57KO"] <- covdata5_reaarangedCPM_sortedavg$ZFP57KO_JB1JF1/(covdata5_reaarangedCPM_sortedavg$ZFP57KO_JB1B6 + covdata5_reaarangedCPM_sortedavg$ZFP57KO_JB1JF1)
head(covdata5_reaarangedCPM_sortedavg)
tail(covdata5_reaarangedCPM_sortedavg,31)
boxplot(tail(covdata5_reaarangedCPM_sortedavg,31)[,c(2:5)]) #Check chrX genes expression
write.table(covdata5_reaarangedCPM_sortedavg, "covdata5_reaarangedCPM_sortedavg.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)

boxplot((covdata5_reaarangedCPM_sortedavg$matRatioWT),
        (covdata5_reaarangedCPM_sortedavg$matRatioZFP57KO), col=(c("darkgreen","darkred")))


#save as covdata5_reaarangedCPM_sortedavg.svg
covdata5_reaarangedCPM_sortedavgbins <- data.frame(cbind(rownames(covdata5_reaarangedCPM_sortedavg), covdata5_reaarangedCPM_sortedavg[,c(5,6)]))
head(covdata5_reaarangedCPM_sortedavgbins)
colnames(covdata5_reaarangedCPM_sortedavgbins) <- c("bins","matRatioWT","matRatioZFP57KO")
head(covdata5_reaarangedCPM_sortedavgbins)
write.table(covdata5_reaarangedCPM_sortedavgbins, 
            "/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb/covdata5_reaarangedCPM_sortedavgbins.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)
#very importat mm10_5mb_%ext.txt created using all possible bins 0-200Mb for all chromosomes (venny intersection), END bins were replced manually to put the correct bins 
#First merge all bins save as merged5mballlbins.txt and replace manually right coordinates end of end points save as merged5mballlbins_replacedright.txt
#sort -k1,1 -k2,2n merged5mballlbins_replacedright.txt > merged5mballlbins_replacedright_sorted.txt

##awk '{print $1"%"$2"%"$3}' merged5mballlbins_replacedright_sorted.txt > mm10_5mb_%ext.txt

mm10_5mb_per <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb/mm10_5Mb_%ext.txt",
                           header=TRUE)
head(mm10_5mb_per)
colnames(mm10_5mb_per) <- "bins"
aggregate_covdata5_per_countedbins = merge(covdata5_reaarangedCPM_sortedavgbins, mm10_5mb_per, by="bins", all.y=TRUE)
head(aggregate_covdata5_per_countedbins)
dim(aggregate_covdata5_per_countedbins)
library(splitstackshape)
covdata5_reaarangedCPM_sortedavgbinsresep <- cSplit(aggregate_covdata5_per_countedbins, "bins", "%")
head(covdata5_reaarangedCPM_sortedavgbinsresep)
covdata5_reaarangedCPM_sortedavgbinsresep <- covdata5_reaarangedCPM_sortedavgbinsresep[,c(3:5,1:2)]
head(covdata5_reaarangedCPM_sortedavgbinsresep)
dim(covdata5_reaarangedCPM_sortedavgbinsresep)
#sort column 1 and 2 same as sort -k1,1 -k2,2n but using R
aggregate_covdata5_per_countedbinsresep_sort <- covdata5_reaarangedCPM_sortedavgbinsresep[order(covdata5_reaarangedCPM_sortedavgbinsresep$bins_1, covdata5_reaarangedCPM_sortedavgbinsresep$bins_2),]
head(aggregate_covdata5_per_countedbinsresep_sort,40)
tail(aggregate_covdata5_per_countedbinsresep_sort)
dim(aggregate_covdata5_per_countedbinsresep_sort)
#Replace NA with -1
aggregate_covdata5_per_countedbinsresep_sort[is.na(aggregate_covdata5_per_countedbinsresep_sort)] <- -1
head(aggregate_covdata5_per_countedbinsresep_sort,110)
write.table(aggregate_covdata5_per_countedbinsresep_sort, 
            "/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb/aggregate_covdata5_per_countedbinsresep_sort", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)
#paste merged5mballlbins_replacedright_sorted.txt in excel add windows 0-n manually and save as mm10_5mb_%ext_bin.txt
mm10_5mb_perbins <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb/mm10_5Mb_%ext_bin.txt",
                               header=TRUE)
#Recombination Map
head(mm10_5mb_perbins,20)
aggregate_covdata5_per_countedbinsresep_sortbin <- data.frame(cbind(mm10_5mb_perbins,aggregate_covdata5_per_countedbinsresep_sort))
head(aggregate_covdata5_per_countedbinsresep_sortbin,20)
write.table(aggregate_covdata5_per_countedbinsresep_sortbin, 
            "/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb/aggregate_covdata5_per_countedbinsresep_sortbin.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)

aggregate_covdata5_per_countedbinsresep_sortbinWT <- aggregate_covdata5_per_countedbinsresep_sortbin[,c(1,5,8)]

TABwt <- xtabs(matRatioWT~., aggregate_covdata5_per_countedbinsresep_sortbinWT)
library(tidyr)
TABwt <- spread(aggregate_covdata5_per_countedbinsresep_sortbinWT, bins_1, matRatioWT, fill=0)
TABwt <- TABwt[,c(1:2,13:20,3:12,21:22)] #rearange column in chr order
head(TABwt)
TABwt[,1]
TABwt <- TABwt[c(1,31,3,14,24:30,32:40,2,4:13,15:23),] #reaarange rows in bins order
rownames(TABwt)=TABwt[,1]
TABwt1 = as.matrix(TABwt[,-1])
colfunc <- colorRampPalette(c("white","yellow","#4169e1","red"))
library(gplots)
library(ggplot2)
library(pheatmap)
svg(filename="TABwt1.5mb.svg", width=10, height=10, pointsize=12)
heatmap.2(TABwt1, Colv = "NA",Rowv ='NA',trace = "none", col = colfunc ,
          lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5),
          density.info=c("none"), dendrogram="none", scale = "none", 
          sepwidth=c(0.001, 0.001), cexRow=1, font=3, cexCol = 0.8, 
          sepcolor="white", margins =c(6,8), srtCol = 45, 
          breaks = c(-1,0,0.31,0.67,1), colsep=1:ncol(TABwt1), rowsep=1:nrow(TABwt1))
dev.off()

aggregate_covdata5_per_countedbinsresep_sortbinZFP57KO <- aggregate_covdata5_per_countedbinsresep_sortbin[,c(1,5,9)]
TABzfp57ko <- xtabs(matRatioZFP57KO~., aggregate_covdata5_per_countedbinsresep_sortbinZFP57KO)
library(tidyr)
TABzfp57ko <- spread(aggregate_covdata5_per_countedbinsresep_sortbinZFP57KO, bins_1, matRatioZFP57KO, fill=0)
TABzfp57ko <- TABzfp57ko[,c(1:2,13:20,3:12,21:22)]
TABzfp57ko <- TABzfp57ko[c(1,31,3,14,24:30,32:40,2,4:13,15:23),]
rownames(TABzfp57ko)=TABzfp57ko[,1]
TABzfp57ko1 = as.matrix(TABzfp57ko[,-1])
colfunc <- colorRampPalette(c("white","yellow","#4169e1","red"))
svg(filename="TABzfp57ko1.5mb.svg", width=10, height=10, pointsize=12)
heatmap.2(TABzfp57ko1, Colv = "NA",Rowv ='NA',trace = "none", col = colfunc ,
          lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5),
          density.info=c("none"), dendrogram="none", scale = "none", 
          sepwidth=c(0.001, 0.001), cexRow=1, font=3, cexCol = 0.8, 
          sepcolor="white", margins =c(6,8), srtCol = 45, 
          breaks = c(-1,0,0.31,0.67,1), colsep=1:ncol(TABzfp57ko1), rowsep=1:nrow(TABzfp57ko1))
dev.off()

head(aggregate_covdata5_per_countedbinsresep_sortbin)
write.table(aggregate_covdata5_per_countedbinsresep_sortbin, "aggregate_covdata5_per_countedbinsresep_sortbin.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE,  col.names = FALSE)

wINDOWS <- data.frame(aggregate_covdata5_per_countedbinsresep_sortbin)
library(tidyverse)
wINDOWS_matRatioWT <-filter(wINDOWS, wINDOWS$matRatioWT > 0.31 & wINDOWS$matRatioWT < 0.67)
head(wINDOWS_matRatioWT)
dim(wINDOWS_matRatioWT)
wINDOWS_matRatioWT <- wINDOWS_matRatioWT[,c(2:4,1,8)]
head(wINDOWS_matRatioWT)
write.table(wINDOWS_matRatioWT, "wINDOWS_matRatioWT.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE,  col.names = FALSE)


wINDOWS_matRatioZFP57KO <-filter(wINDOWS, wINDOWS$matRatioZFP57KO > 0.31 & wINDOWS$matRatioZFP57KO < 0.67)
head(wINDOWS_matRatioZFP57KO)
dim(wINDOWS_matRatioZFP57KO)
wINDOWS_matRatioZFP57KO <- wINDOWS_matRatioZFP57KO[,c(2:4,1,9)]
head(wINDOWS_matRatioZFP57KO)
write.table(wINDOWS_matRatioZFP57KO, "wINDOWS_matRatioZFP57KO.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE, col.names = FALSE)

NonwINDOWS <- data.frame(aggregate_covdata5_per_countedbinsresep_sortbin)
library(tidyverse)
NonwINDOWS_matRatioWT <-filter(NonwINDOWS, NonwINDOWS$matRatioWT < 0.31 | NonwINDOWS$matRatioWT > 0.67 |  NonwINDOWS$matRatioWT < 0)
head(NonwINDOWS_matRatioWT)
dim(NonwINDOWS_matRatioWT)
NonwINDOWS_matRatioWT <- NonwINDOWS_matRatioWT[,c(2:4,1,8)]
head(NonwINDOWS_matRatioWT)
write.table(NonwINDOWS_matRatioWT, "NonwINDOWS_matRatioWT.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE,  col.names = FALSE)


NonwINDOWS_matRatioZFP57KO <-filter(NonwINDOWS, NonwINDOWS$matRatioZFP57KO < 0.31 | NonwINDOWS$matRatioZFP57KO > 0.67 |  NonwINDOWS$matRatioWT < 0)
head(NonwINDOWS_matRatioZFP57KO)
dim(NonwINDOWS_matRatioZFP57KO)
NonwINDOWS_matRatioZFP57KO <- NonwINDOWS_matRatioZFP57KO[,c(2:4,1,9)]
head(NonwINDOWS_matRatioZFP57KO)
write.table(NonwINDOWS_matRatioZFP57KO, "NonwINDOWS_matRatioZFP57KO.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE, col.names = FALSE)
#WT with ICR genes
sort -k1,1 -k2,2n wINDOWS_matRatioWT.txt > wINDOWS_matRatioWTsort.txt
bedtools intersect -wa -wb -a wINDOWS_matRatioWTsort.txt -b /home/ankitv/ref_av/mm10/mm10_ICR.bed 
sort -k1,1 -k2,2n NonwINDOWS_matRatioWT.txt > NonwINDOWS_matRatioWTsort.txt
bedtools intersect -wa -wb -a NonwINDOWS_matRatioWTsort.txt -b /home/ankitv/ref_av/mm10/mm10_ICR.bed 

#ZFP57KO with ICR genes
sort -k1,1 -k2,2n wINDOWS_matRatioZFP57KO.txt > wINDOWS_matRatioZFP57KOsort.txt
bedtools intersect -wa -wb -a wINDOWS_matRatioZFP57KOsort.txt -b /home/ankitv/ref_av/mm10/mm10_ICR.bed 
sort -k1,1 -k2,2n NonwINDOWS_matRatioZFP57KO.txt > NonwINDOWS_matRatioZFP57KOsort.txt
bedtools intersect -wa -wb -a NonwINDOWS_matRatioZFP57KOsort.txt -b /home/ankitv/ref_av/mm10/mm10_ICR.bed 

#WT with Imprinted genes
#sort -k1,1 -k2,2n wINDOWS_matRatioWT.txt > wINDOWS_matRatioWTsort.txt
##bedtools intersect -wa -wb -a wINDOWS_matRatioWTsort.txt -b imprinted_gene_cords.bed
#sort -k1,1 -k2,2n NonwINDOWS_matRatioWT.txt > NonwINDOWS_matRatioWTsort.txt
##bedtools intersect -wa -wb -a NonwINDOWS_matRatioWTsort.txt -b imprinted_gene_cords.bed

#ZFP57KO with Imprinted genes
#sort -k1,1 -k2,2n wINDOWS_matRatioZFP57KO.txt > wINDOWS_matRatioZFP57KOsort.txt
##bedtools intersect -wa -wb -a wINDOWS_matRatioZFP57KOsort.txt -b imprinted_gene_cords.bed
#sort -k1,1 -k2,2n NonwINDOWS_matRatioZFP57KO.txt > NonwINDOWS_matRatioZFP57KOsort.txt
##bedtools intersect -wa -wb -a NonwINDOWS_matRatioZFP57KOsort.txt -b imprinted_gene_cords.bed

#Filter non-windows/biased regions : For detailed description see Oct2019 in documents folder
cat NonwINDOWS_matRatioWTsort.txt NonwINDOWS_matRatioZFP57KOsort.txt > combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort.txt
awk '{if($5 >= 0) print}' combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort.txt > combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort_coverage.txt 
awk '{if($5 < 0) print}' combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort.txt > combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort_NOcoverage.txt
awk '{print $1"%"$2"%"$3}' combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort_coverage.txt | sort -k1,1 -u | awk -F'%' '{print $1"\t"$2"\t"$3}' | sort -k1,1 -k2,2n > combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort_coverage_toremove.txt



#Import Bulk data and calculate size factor (Deseq2)
##################################################################################### DESeq2 ###############################################################################################################
library(DESeq2)
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2")
Bulkcountdata <- read.table("JB1_time12_star-featureCounts_mm10_Igf2.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
Bulkcountdata <- Bulkcountdata[ ,6:ncol(Bulkcountdata)]
head(Bulkcountdata)
# Remove .bam or .sam from filenames
colnames(Bulkcountdata) <- gsub("\\.[sb]am$", "", colnames(Bulkcountdata))
dim(Bulkcountdata)
head(Bulkcountdata)
Bulkcountdata = as.matrix(Bulkcountdata)
head(Bulkcountdata)
colSums(Bulkcountdata) #7663520 10225378  9545119 10337082 > use this for CPM normalization
Bulkcoldata <- read.table("Bulkcoldata.txt" , header = TRUE, stringsAsFactors = FALSE)
rownames(Bulkcoldata)
rownames(Bulkcoldata) <- colnames(Bulkcountdata)
head(Bulkcoldata)
Bulkcoldata <- Bulkcoldata[,c("condition","replicate")]
Bulkcoldata$condition <- factor(Bulkcoldata$condition)
Bulkcoldata$replicate <- factor(Bulkcoldata$replicate)
all(rownames(Bulkcoldata) == colnames(Bulkcountdata)) #should print TRUE

Bulkdds <- DESeqDataSetFromMatrix(countData=Bulkcountdata, colData = Bulkcoldata, design = ~ condition)
Bulkdds
Bulkkeep <- rowSums(counts(Bulkdds)) >= 9
Bulkdds <- Bulkdds[Bulkkeep,]
Bulkddscounts <- counts(Bulkdds, normalized=FALSE)
head(Bulkddscounts)
dim(Bulkddscounts)
colSums(Bulkddscounts)

#NOTE: For CPM calculations use original bulkcounts sum (non filtered  for rowsum >9) but for deseq2 normalization size factor should be calculated after filt and use the same

#View filtered count matrix: View(counts(Bulkdds))
#Normalization is the part of DESeq command: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#Normalized separately:
BulkddsNorm <- estimateSizeFactors(Bulkdds)
sizeFactors(BulkddsNorm) #0.8706398  1.1225115  1.0022048 1.0605409


#Note: it should be noted that gtf used for bulk were cleaned for scaffolds and then used for allele-specific 
#..as the chromosome name need to be changed from  chrN to N, scaffolds doesnot have chr only alphabets
# Therefore features in bulk are 55604 and allele sp has 54750
#Igf2 matrix removed ENSMUSG00000070140.1  ENSMUSG00000115302.1 from  gtf. These genes shared a large coordinate range with Igf2 and therefore impair read assignment to Igf2, subtracted gtf file is gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschrlike-.gtf
#fgrep -f igf2likegenes.txt ./../../N-masked-JF1-GRCm38-M20-overlapped/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschr.gtf -v > gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschrlike-.gtf
#Note sorted by readname
#/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -p -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-GRCm38-M20-overlapped/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschrlike-.gtf -T 12 -o ./Igf2/Alsp_JB1_time12_star-featureCounts_GRCm38.mm10_Igf2.txt JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam


#Featurecount beforesplit Igf2
#/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -p -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-GRCm38-M20-overlapped/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschrlike-.gtf -T 12 -o Alsp_JB1_time12_star-featureCounts_GRCm38.mm10.BEFORESPLIT_Igf2.txt JB1_WT_Rep1_Aligned.sortedByCoord.out.bam JB1_WT_Rep2_Aligned.sortedByCoord.out.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bam

#Igf2 matrix removed ENSMUSG00000070140.1  ENSMUSG00000115302.1 from  gtf. These genes shared a large coordinate range with Igf2 and therefore impair read assignment to Igf2, subtracted gtf file is gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minusigf2like.gtf

#fgrep -f igf2likegenes.txt /home/ankitv/ref_av/gencodes/gencode_M20/prep/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf -v > gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minusigf2like.gtf
#Import Before split data and calculate library size and size factor (DESeq2)

############################################# Featurecount ##############################################
#/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -p -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-Grcm38_overlapped/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschr.gtf -T 12 -o Alsp_JB1_time12_star-featureCounts_GRCm38.txt JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam
#/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -p -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-Grcm38_overlapped/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschrlike-.gtf -T 12 -o ./Igf2/Alsp_JB1_time12_star-featureCounts_GRCm38.mm10_Igf2.txt JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam

################################################################# DESeq2 ######################################################################
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb")
gcountdata <- read.table("Alsp_JB1_time12_star-featureCounts_GRCm38.mm10_Igf2.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
gcountdata <- gcountdata[ ,6:ncol(gcountdata)]
head(gcountdata)
# Remove .bam or .sam from filenames
colnames(gcountdata) <- gsub("\\.[sb]am$", "", colnames(gcountdata))
colnames(gcountdata) <-  c("JB1_WT_Rep1.B6","JB1_WT_Rep2.B6","JB1_WT_Rep1.JF1","JB1_WT_Rep2.JF1","JB1_ZFP57_KO_Rep1.B6","JB1_ZFP57_KO_Rep2.B6","JB1_ZFP57_KO_Rep1.JF1","JB1_ZFP57_KO_Rep2.JF1")
head(gcountdata)
dim(gcountdata)
head(gcountdata)
library(DESeq2)
gcountdata = as.matrix(gcountdata)
head(gcountdata)
boxplot(gcountdata, ylim=c(0,1))
plotPCA(as.matrix(gcountdata))  #library(EDASeq)
plotPCA(as.matrix(gcountdata), labels=F, col =  c("blue","blue","darkred","darkred","darkgreen","darkgreen","orange","orange"))

ggsave("PCA_tgcountdatafilt_scaleT.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

gcoldata <- read.table("gcoldata.txt" , header = TRUE, stringsAsFactors = FALSE)
rownames(gcoldata)
gcoldata[,1]
rownames(gcoldata)=gcoldata[,1]
rownames(gcoldata)
colnames(gcoldata)
gcoldata = gcoldata[,-1]
gcoldata <- data.frame(gcoldata)
head(gcoldata)
gcoldata <- gcoldata[,c("gcondition","greplicate","gallele")]
rownames(gcoldata) <- colnames(gcountdata)
gcoldata$gcondition <- factor(gcoldata$gcondition)
all(rownames(gcoldata) == colnames(gcountdata)) #should print TRUE
gdds <- DESeqDataSetFromMatrix(countData =gcountdata, colData = gcoldata, design = ~ gcondition)
gdds
#featureData <- data.frame(gene=rownames(gcountdata))
keep <- rowSums(counts(gdds)) >= 10
gdds <- gdds[keep,]
#View filtered count matrix: View(counts(gdds))
#Apply first method Size factor based of DESeq2
gddsSF <- gdds
gddsSF <- estimateSizeFactors(gddsSF)
sizeFactors(gddsSF)
#Note normalization =FALSE will not divide counts by normalization factors
gddsSFcounts <- counts(gddsSF, normalized=FALSE)
head(gddsSFcounts)
boxplot(gddsSFcounts, ylim=c(0,100))
plotPCA(as.matrix(gddsSFcounts))
plotPCA(as.matrix(gddsSFcounts), labels=F, col =  c("blue","blue","darkred","darkred","darkgreen","darkgreen","orange","orange"))
tgddsSFcountsfilt = t(gddsSFcounts)
dim(tgddsSFcountsfilt)
rownames(tgddsSFcountsfilt)
tgddsSFcountsfilt = data.frame(tgddsSFcountsfilt)
tgddsSFcountsfilt["Color"] <-  c("WTB6_1","WTB6_2","WTJF1_1","WTJF1_2","ZFP57KOB6_1","ZFP57KOB6_2","ZFP57KOJF1_1","ZFP57KOJF1_2")
dim(tgddsSFcountsfilt)
dfx <-tgddsSFcountsfilt[c(1:15673)]
PC<-prcomp(dfx, scale. = T)
PCi<-data.frame(PC$x,Color=tgddsSFcountsfilt$Color)
percentage <- round(PC$sdev^2 / sum(PC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentage <- paste( colnames(PCi), "(", paste( as.character(percentage), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

p12<-ggplot(PCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentage[1]) + ylab(percentage[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("blue","blue","darkred","darkred","darkgreen","darkgreen","orange","orange"))+
  scale_shape_manual(values=c(19,17,19,17,19,17,19,17))
p12 <- p12+theme_bw()
p12
#save as p12.svg
#Normalization is the part of DESeq command: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#But here I set size factor as per bulk  counts
#sizeFactors(BulkddsNorm) #0.8706398  1.1225115  1.0022048 1.0605409
gddsSFNorm <- gdds
sizeFactors(gddsSFNorm) <- c(0.8706398,1.1225115,0.8706398,1.1225115,1.0022048,1.0605409,1.0022048,1.0605409)
sizeFactors(gddsSFNorm)
#Normalize by DESeq2 size factor, set size factoe using Bulkdds 
#Note normalization =TRUE divide counts by the user set normalization factors
#Normalize allele specfic counts export normalized counts
gddsSFNormcounts <- counts(gddsSFNorm, normalized=TRUE)
head(gddsSFNormcounts)
write.table(gddsSFNormcounts, "gddsSFNormcounts.txt", sep="\t", quote=F, col.names=NA)
boxplot(gddsSFNormcounts, ylim=c(0,100))
plotPCA(as.matrix(gddsSFNormcounts))
plotPCA(as.matrix(gddsSFNormcounts), labels=F, col =  c("blue","blue","darkred","darkred","darkgreen","darkgreen","orange","orange"))
tgddsSFNormcountsfilt = t(gddsSFNormcounts)
dim(tgddsSFNormcountsfilt)
rownames(tgddsSFNormcountsfilt)
tgddsSFNormcountsfilt = data.frame(tgddsSFNormcountsfilt)
tgddsSFNormcountsfilt["Color"] <-  c("WTB6_1","WTB6_2","WTJF1_1","WTJF1_2","ZFP57KOB6_1","ZFP57KOB6_2","ZFP57KOJF1_1","ZFP57KOJF1_2")
dim(tgddsSFNormcountsfilt)
dfx <-tgddsSFNormcountsfilt[c(1:15673)]
PC<-prcomp(dfx, scale. = T)
PCi<-data.frame(PC$x,Color=tgddsSFNormcountsfilt$Color)
percentage <- round(PC$sdev^2 / sum(PC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentage <- paste( colnames(PCi), "(", paste( as.character(percentage), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

p13<-ggplot(PCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentage[1]) + ylab(percentage[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("blue","blue","darkred","darkred","darkgreen","darkgreen","orange","orange"))+
  scale_shape_manual(values=c(19,17,19,17,19,17,19,17))
p13 <- p13+theme_bw()
p13
ggsave("PCA_tgddsSFNormcountsfilt_scaleT.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T

#z-score calc
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb")
gddsSFNormcounts <- read.table("gddsSFNormcounts.txt" , header = TRUE, stringsAsFactors = FALSE)
head(gddsSFNormcounts)
dim(gddsSFNormcounts)
gddsSFNormcounts <- data.frame(gddsSFNormcounts[order(rownames(gddsSFNormcounts)),])
head(gddsSFNormcounts)
#Transform and scale
z_TgddsSFNormcounts = scale(t(gddsSFNormcounts), center = TRUE, scale = TRUE)
z_gddsSFNormcounts <- t(z_TgddsSFNormcounts)
head(z_gddsSFNormcounts)
write.table(z_gddsSFNormcounts, "z_gddsSFNormcounts.txt", sep="\t", quote=F, col.names=NA)

#z-score calc on averaged normalized counts
head(gddsSFNormcounts)
dim(gddsSFNormcounts)
gddsSFNormcounts <- data.frame(gddsSFNormcounts)
#gddsSFNormcountsavg <- data.frame(cbind(rownames(gddsSFNormcounts),(gddsSFNormcounts$JB1_WT_Rep1.B6+gddsSFNormcounts$JB1_WT_Rep2.B6)/2, (gddsSFNormcounts$JB1_WT_Rep1.JF1 + gddsSFNormcounts$JB1_WT_Rep2.JF1)/2, (gddsSFNormcounts$JB1_ZFP57_KO_Rep1.B6 + gddsSFNormcounts$JB1_ZFP57_KO_Rep2.B6)/2, (gddsSFNormcounts$JB1_ZFP57_KO_Rep1.JF1 + gddsSFNormcounts$JB1_ZFP57_KO_Rep2.JF1)/2))
gddsSFNormcountsavg <- data.frame(cbind(data.frame(rownames(gddsSFNormcounts)),
                                        (gddsSFNormcounts$JB1_WT_Rep1.B6 + gddsSFNormcounts$JB1_WT_Rep2.B6)/2,
                                        (gddsSFNormcounts$JB1_WT_Rep1.JF1 + gddsSFNormcounts$JB1_WT_Rep2.JF1)/2, 
                                        (gddsSFNormcounts$JB1_ZFP57_KO_Rep1.B6 + gddsSFNormcounts$JB1_ZFP57_KO_Rep2.B6)/2, 
                                        (gddsSFNormcounts$JB1_ZFP57_KO_Rep1.JF1 + gddsSFNormcounts$JB1_ZFP57_KO_Rep2.JF1)/2))
head(gddsSFNormcountsavg)
rownames(gddsSFNormcountsavg) <- gddsSFNormcountsavg[,1]
gddsSFNormcountsavg <- gddsSFNormcountsavg[,-1]
head(gddsSFNormcountsavg)
dim(gddsSFNormcountsavg)
colnames(gddsSFNormcountsavg) <- c("JB1_WTB6","JB1_WTJF1","JB1_ZFP57KOB6","JB1_ZFP57KOJF1")
rownames(gddsSFNormcountsavg)
head(gddsSFNormcountsavg)
#Transpose and scale
z_TgddsSFNormcountsavg= scale(t(gddsSFNormcountsavg), center = TRUE, scale = TRUE)
z_gddsSFNormcountsavg <- t(z_TgddsSFNormcountsavg)
head(z_gddsSFNormcountsavg)
dim(z_gddsSFNormcountsavg)
write.table(z_gddsSFNormcountsavg, "z_gddsSFNormcountsavg.txt", sep="\t", quote=F, col.names=NA)


#gddsa
gcoldata$gallele <- factor(gcoldata$gallele)
gddsa <- DESeqDataSetFromMatrix(countData =gcountdata, colData = gcoldata, design = ~ gallele)
gddsa
#featureData <- data.frame(gene=rownames(gcountdata))
keep <- rowSums(counts(gddsa)) >= 10
gddsa <- gddsa[keep,]

#View filtered count matrix: View(counts(gdds))
#Apply first method Size factor based of DESeq2
gddsaSF <- gddsa
gddsaSF <- estimateSizeFactors(gddsaSF)
sizeFactors(gddsaSF)
#Note normalization =FALSE will not divide counts by normalization factors
gddsaSFcounts <- counts(gddsaSF, normalized=FALSE)
head(gddsaSFcounts)
gddsaSFcounts <- data.frame(gddsaSFcounts[order(rownames(gddsaSFcounts)),])
head(gddsaSFcounts)
#Normalization is the part of DESeq command: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html

gddsaSFNorm <- gddsa
#Normalize by DESeq2 size factor, set size factoe using Bulkdds 
sizeFactors(gddsaSFNorm) <- c(0.8706398,1.1225115,0.8706398,1.1225115,1.0022048,1.0605409,1.0022048,1.0605409)
sizeFactors(gddsaSFNorm)
#Note normalization =TRUE divide counts by the user set normalization factors
#Normalize allele specfic counts: export normalized counts
gddsaSFNormcounts <- counts(gddsaSFNorm, normalized=TRUE)
head(gddsaSFNormcounts)
write.table(gddsaSFNormcounts, "gddsaSFNormcounts.txt", sep="\t", quote=F)
boxplot(gddsaSFNormcounts, ylim=c(0,100))

tgddsaSFNormcountsfilt = t(gddsaSFNormcounts)
dim(tgddsaSFNormcountsfilt)
rownames(tgddsaSFNormcountsfilt)
tgddsaSFNormcountsfilt = data.frame(tgddsaSFNormcountsfilt)
tgddsaSFNormcountsfilt["Color"] <-  c("WTB6_1","WTB6_2","WTJF1_1","WTJF1_2","ZFP57KOB6_1","ZFP57KOB6_2","ZFP57KOJF1_1","ZFP57KOJF1_2")
dim(tgddsaSFNormcountsfilt)
dfx <-tgddsaSFNormcountsfilt[c(1:15673)]
PC<-prcomp(dfx, scale. = T)
PCi<-data.frame(PC$x,Color=tgddsaSFNormcountsfilt$Color)
percentage <- round(PC$sdev^2 / sum(PC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentage <- paste( colnames(PCi), "(", paste( as.character(percentage), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

p15<-ggplot(PCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentage[1]) + ylab(percentage[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("blue","blue","darkred","darkred","darkgreen","darkgreen","orange","orange"))+
  scale_shape_manual(values=c(19,17,19,17,19,17,19,17))
p15 <- p15+theme_bw()
p15

#save as PCA_gddsaSFNormcounts.svg
#z-score calc: gddsa
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb")
gddsaSFNormcounts <- read.table("gddsaSFNormcounts.txt" , header = TRUE, stringsAsFactors = FALSE)
head(gddsaSFNormcounts)
dim(gddsaSFNormcounts)
gddsaSFNormcounts <- data.frame(gddsaSFNormcounts[order(rownames(gddsaSFNormcounts)),])
head(gddsaSFNormcounts)
#Transform and scale
z_TgddsaSFNormcounts= scale(t(gddsaSFNormcounts), center = TRUE, scale = TRUE)
z_gddsaSFNormcounts <- t(z_TgddsaSFNormcounts)
head(z_gddsaSFNormcounts)
write.table(z_gddsaSFNormcounts, "z_gddsaSFNormcounts.txt", sep="\t", quote=F, col.names=NA)

#z-score calc on averaged normalized counts
head(gddsaSFNormcounts)
dim(gddsaSFNormcounts)
gddsaSFNormcounts <- data.frame(gddsaSFNormcounts)
#gddsaSFNormcountsavg <- data.frame(cbind(rownames(gddsaSFNormcounts),(gddsaSFNormcounts$JB1_WT_Rep1.B6+gddsaSFNormcounts$JB1_WT_Rep2.B6)/2, (gddsaSFNormcounts$JB1_WT_Rep1.JF1 + gddsaSFNormcounts$JB1_WT_Rep2.JF1)/2, (gddsaSFNormcounts$JB1_ZFP57_KO_Rep1.B6 + gddsaSFNormcounts$JB1_ZFP57_KO_Rep2.B6)/2, (gddsaSFNormcounts$JB1_ZFP57_KO_Rep1.JF1 + gddsaSFNormcounts$JB1_ZFP57_KO_Rep2.JF1)/2))
gddsaSFNormcountsavg <- data.frame(cbind(data.frame(rownames(gddsaSFNormcounts)),
                                         (gddsaSFNormcounts$JB1_WT_Rep1.B6 + gddsaSFNormcounts$JB1_WT_Rep2.B6)/2,
                                         (gddsaSFNormcounts$JB1_WT_Rep1.JF1 + gddsaSFNormcounts$JB1_WT_Rep2.JF1)/2, 
                                         (gddsaSFNormcounts$JB1_ZFP57_KO_Rep1.B6 + gddsaSFNormcounts$JB1_ZFP57_KO_Rep2.B6)/2, 
                                         (gddsaSFNormcounts$JB1_ZFP57_KO_Rep1.JF1 + gddsaSFNormcounts$JB1_ZFP57_KO_Rep2.JF1)/2))
head(gddsaSFNormcountsavg)
rownames(gddsaSFNormcountsavg) <- gddsaSFNormcountsavg[,1]
gddsaSFNormcountsavg <- gddsaSFNormcountsavg[,-1]
head(gddsaSFNormcountsavg)
dim(gddsaSFNormcountsavg)
colnames(gddsaSFNormcountsavg) <- c("JB1_WTB6","JB1_WTJF1","JB1_ZFP57KOB6","JB1_ZFP57KOJF1")
rownames(gddsaSFNormcountsavg)
head(gddsaSFNormcountsavg)
gddsaSFNormcountsavg.id <- data.frame(rownames(gddsaSFNormcountsavg))
colnames(gddsaSFNormcountsavg.id) <- "id"
head(gddsaSFNormcountsavg.id)

#Transform and scale
z_TgddsaSFNormcountsavg= scale(t(gddsaSFNormcountsavg), center = TRUE, scale = TRUE)
z_gddsaSFNormcountsavg <- t(z_TgddsaSFNormcountsavg)
head(z_gddsaSFNormcountsavg)
dim(z_gddsaSFNormcountsavg)
write.table(z_gddsaSFNormcountsavg, "z_gddsaSFNormcountsavg.txt", sep="\t", quote=F, col.names=NA)
z_gddsaSFNormcountsavg.id <- data.frame(cbind(gddsaSFNormcountsavg.id, z_gddsaSFNormcountsavg))
head(z_gddsaSFNormcountsavg.id)
tail(z_gddsaSFNormcountsavg.id)
write.table(z_gddsaSFNormcountsavg.id, "z_gddsaSFNormcountsavg.id.txt", sep="\t", quote=F, col.names=T, row.names = F)

chr.pos =  read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt",header=FALSE)
head(chr.pos)
colnames(chr.pos) <- c("id", "Genes", "chr", "start", "end")
head(chr.pos)
chr.pos.2 = merge(gddsaSFNormcountsavg.id, chr.pos, by="id", all.x=TRUE)
head(chr.pos.2)
dim(chr.pos.2)
gddsaSFNormcountsavg_chr_gene <- data.frame(cbind(chr.pos.2, gddsaSFNormcountsavg))
head(gddsaSFNormcountsavg_chr_gene)
gddsaSFNormcountsavg_chr_gene1 <- gddsaSFNormcountsavg_chr_gene[,c(3:5,2,1,6:9)]
head(gddsaSFNormcountsavg_chr_gene1)
write.table(gddsaSFNormcountsavg_chr_gene1, "gddsaSFNormcountsavg_chr_gene1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)

chr.pos.3 = merge(z_gddsaSFNormcountsavg.id, chr.pos, by="id", all.x=TRUE)
head(chr.pos.3)
dim(chr.pos.3)
z_gddsaSFNormcountsavg_chr_gene <- data.frame(cbind(chr.pos.3, z_gddsaSFNormcountsavg))
head(z_gddsaSFNormcountsavg_chr_gene)
z_gddsaSFNormcountsavg_chr_gene1 <- z_gddsaSFNormcountsavg_chr_gene[,c(7:9,6,1,2:5)]
head(z_gddsaSFNormcountsavg_chr_gene1)
tail(z_gddsaSFNormcountsavg_chr_gene1)
dim(z_gddsaSFNormcountsavg_chr_gene1)
write.table(z_gddsaSFNormcountsavg_chr_gene1, "z_gddsaSFNormcountsavg_chr_gene1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)


#Individuals replicates separate process
head(gddsaSFNormcounts)
gddsaSFNormcounts["id"] <- data.frame(rownames(gddsaSFNormcounts))
head(gddsaSFNormcounts)
chr.pos.4 = merge(gddsaSFNormcounts, chr.pos, by="id", all.x=TRUE)
head(chr.pos.4)
dim(chr.pos.4)
gddsaSFNormcounts_chr_gene <- chr.pos.4
head(gddsaSFNormcounts_chr_gene)
gddsaSFNormcounts_chr_gene1 <- gddsaSFNormcounts_chr_gene[,c(11:13,10,1,2:9)]
head(gddsaSFNormcounts_chr_gene1)
write.table(gddsaSFNormcounts_chr_gene1, "gddsaSFNormcounts_chr_gene1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)

gddsaSFNormcounts_chr_gene1 <- read.table("gddsaSFNormcounts_chr_gene1.txt", header = T, stringsAsFactors = F)
head(gddsaSFNormcounts_chr_gene1,1)
#All  genes mat ratio
gddsaSFNormcounts_chr_gene1_ratio <- gddsaSFNormcounts_chr_gene1
gddsaSFNormcounts_chr_gene1_ratio["WT_allele_ratio_1"] <- gddsaSFNormcounts_chr_gene1_ratio$JB1_WT_Rep1.JF1 / (gddsaSFNormcounts_chr_gene1_ratio$JB1_WT_Rep1.B6 + gddsaSFNormcounts_chr_gene1_ratio$JB1_WT_Rep1.JF1)
gddsaSFNormcounts_chr_gene1_ratio["WT_allele_ratio_2"] <- gddsaSFNormcounts_chr_gene1_ratio$JB1_WT_Rep2.JF1 / (gddsaSFNormcounts_chr_gene1_ratio$JB1_WT_Rep2.B6 + gddsaSFNormcounts_chr_gene1_ratio$JB1_WT_Rep2.JF1)
gddsaSFNormcounts_chr_gene1_ratio["WT_allele_ratio"] <- (gddsaSFNormcounts_chr_gene1_ratio$WT_allele_ratio_1 + gddsaSFNormcounts_chr_gene1_ratio$WT_allele_ratio_2)/2
gddsaSFNormcounts_chr_gene1_ratio["ZFP57KO_allele_ratio_1"] <- gddsaSFNormcounts_chr_gene1_ratio$JB1_ZFP57_KO_Rep1.JF1 / (gddsaSFNormcounts_chr_gene1_ratio$JB1_ZFP57_KO_Rep1.B6 + gddsaSFNormcounts_chr_gene1_ratio$JB1_ZFP57_KO_Rep1.JF1)
gddsaSFNormcounts_chr_gene1_ratio["ZFP57KO_allele_ratio_2"] <- gddsaSFNormcounts_chr_gene1_ratio$JB1_ZFP57_KO_Rep2.JF1 / (gddsaSFNormcounts_chr_gene1_ratio$JB1_ZFP57_KO_Rep2.B6 + gddsaSFNormcounts_chr_gene1_ratio$JB1_ZFP57_KO_Rep2.JF1)
gddsaSFNormcounts_chr_gene1_ratio["ZFP57KO_allele_ratio"] <- (gddsaSFNormcounts_chr_gene1_ratio$ZFP57KO_allele_ratio_1 + gddsaSFNormcounts_chr_gene1_ratio$ZFP57KO_allele_ratio_2)/2
head(gddsaSFNormcounts_chr_gene1_ratio)
write.table(gddsaSFNormcounts_chr_gene1_ratio, "gddsaSFNormcounts_chr_gene1_ratio.txt", sep = "\t", quote = F, append = F, row.names = F)

#Get imprinted genes
#This file imprinted_gene_name.dups.txt was produced from my search and bouschet imprinted genes and duplicates were removed
sort -k1,1 imprinted_gene_name.dups.txt -u > imprinted_gene_name.txt #188
/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt
#I tried assigning name manually but very less coordinates were assigned: fgrep -f imprinted_gene_name.txt /home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt -w | wc -l#158
#So I rechecked genes symbols of all genes which were not found using online MGI and assign IDs of MGI or lifover from paper 1 coordinate
#Added everything  in a new file and save as imprinted_gene_rev.txt
#Sort the file chr wise
grep chr imprinted_gene_rev.txt | sort -k1,1 -k2,2n > imprinted_gene_rev.sort.txt #169 genes
#Not used here library(biomaRt)

#ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
#mouse_gene_ids  <- read.table("imprinted_gene_name.txt")
#allimprinted <- getBM(attributes=c('chromosome_name','start_position','end_position','ensembl_gene_id','external_gene_name'),
#                      filters = 'mgi_symbol',
#                      values = mouse_gene_ids,
#                      mart = ensembl)
#allimprinted
#write.table(allimprinted, "allimprinted.txt",sep="\t", quote = F, append = F, row.names = F)  #Ensembl ENSID.version(ENSMUSG00000010751.15)  is not required because version will change 
#awk '{print "chr"$1"\t"$2"\t"$3"\t"$4"\t"$5}' allimprinted.txt > allimprinted_chr.txt

awk '{print $4"\t"$5}' imprinted_gene_rev.sort.txt | grep ENS > imprinted_gene_name_ens.txt
awk '{print $1}' imprinted_gene_name_ens.txt > imprinted_gene_ens.txt  #159 genes 10 were not presemnt in ensembl
#Use fgrep, dont use (-w). This command will give genename and ensembl.version


#Extract imprinted genes data
#Selected genes list was given by Basilia (36) and coordinates were obtained by fgrep -f with gtf gene file. 
#Some genes like Inpp5f and Ddc were filtere out because they were not found to be imprinted leaving 34 genes 
#I removed it before itself as they can affect the z-score scaling so I kept only 34 genes

selectedgenes.cords <- read.table("selected_imprinted_genes_refilt.cords.txt", header = F)
head(selectedgenes.cords)
selectedgenes.cords<- selectedgenes.cords[,c(1:3,5:6)]
colnames(selectedgenes.cords) <- c("chr", "start", "end", "id", "Gene")
z_gddsaSFNormcountsavg.selected_genes = merge(z_gddsaSFNormcountsavg.id, selectedgenes.cords, by="id", all.x=FALSE)
head(z_gddsaSFNormcountsavg.selected_genes)
dim(z_gddsaSFNormcountsavg.selected_genes)
library(ggplot2)
library(gplots)
#setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb")
#data <- read.table("all_DEG_rearranged_filt_ImpGenes.heatmap.txt",header=TRUE)
z_gddsaSFNormcountsavg.selected_genes_sort <- z_gddsaSFNormcountsavg.selected_genes[order(z_gddsaSFNormcountsavg.selected_genes$chr, z_gddsaSFNormcountsavg.selected_genes$start),]
head(z_gddsaSFNormcountsavg.selected_genes_sort)
dim(z_gddsaSFNormcountsavg.selected_genes_sort)
z_gddsaSFNormcountsavg.selected_genes_sort <- z_gddsaSFNormcountsavg.selected_genes_sort[c(1,13:34, 2:12),]
head(z_gddsaSFNormcountsavg.selected_genes_sort)
z_gddsaSFNormcountsavg.selected_genes_sort_heat <- z_gddsaSFNormcountsavg.selected_genes_sort[,c(9,2:5)]
rownames(z_gddsaSFNormcountsavg.selected_genes_sort_heat)
rownames(z_gddsaSFNormcountsavg.selected_genes_sort_heat)=z_gddsaSFNormcountsavg.selected_genes_sort_heat[,1]
rownames(z_gddsaSFNormcountsavg.selected_genes_sort_heat)
colnames(z_gddsaSFNormcountsavg.selected_genes_sort_heat)
head(z_gddsaSFNormcountsavg.selected_genes_sort_heat)
z_gddsaSFNormcountsavg.selected_genes_sort_heat1 = z_gddsaSFNormcountsavg.selected_genes_sort_heat[,-1]
z_gddsaSFNormcountsavg.selected_genes_sort_heat1 = as.matrix(z_gddsaSFNormcountsavg.selected_genes_sort_heat1)
head(z_gddsaSFNormcountsavg.selected_genes_sort_heat1)
dim(z_gddsaSFNormcountsavg.selected_genes_sort_heat1)
write.table(z_gddsaSFNormcountsavg.selected_genes_sort_heat1, "z_gddsaSFNormcountsavg.selected_genes_sort_heat1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)

#colfunc <- colorRampPalette( colors = brewer.pal(9,"PRGn") )
colfunc <- colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","white","#FDDBC7","#F4A582","#D6604D","#B2182B"))

#heatmap.2(data1, Colv = "NA",trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=1, font=3, cexCol = 1, margins =c(3,8), srtCol = 45, breaks = seq(-5,5, length.out = 100))
svg(filename="z_gddsaSFNormcountsavg.selected_genes_sort_heat1.svg", width=8, height=18, pointsize=12)
heatmap.2(z_gddsaSFNormcountsavg.selected_genes_sort_heat1, Colv = "NA",Rowv ='NA',trace = "none", col = colfunc ,
          lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5),
          density.info=c("none"), dendrogram="none", scale = "none", 
          sepwidth=c(0.001, 0.001), cexRow=2, font=3, cexCol = 0.8, 
          sepcolor="black", margins =c(6,8), srtCol = 45, 
          breaks = seq(-1,1, length.out = 100), colsep=1:ncol(z_gddsaSFNormcountsavg.selected_genes_sort_heat1),
          rowsep=1:nrow(z_gddsaSFNormcountsavg.selected_genes_sort_heat1))
dev.off()

library(pheatmap)
library(RColorBrewer)
breaksList = seq(-1, 1)
#pheatmap(data1,treeheight_row = 0, cluster_cols=F, cluster_rows=F, treeheight_col = 0, gaps_col =NULL, gaps_row = NULL, border_color = "black", breaks = breaksList,color= colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)))
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_gddsaSFNormcountsavg.selected_genes_sort_heat1,
         color = colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","white","#FDDBC7","#F4A582","#D6604D","#B2182B"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=F, 
         cluster_rows=F,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12, 
         filename = "z_gddsaSFNormcountsavg.selected_genes_sort_heat1.png")


dev.off()
#Now filter the matrix to remove low expressed genes, so cpm counts >=rowsum 10 are removed
#write.table(results, "results.txt", sep="\t", quote=F, col.names=NA)
gddsaSFNormcountsavg.id.avg <- data.frame(cbind(gddsaSFNormcountsavg.id, gddsaSFNormcountsavg))
head(gddsaSFNormcountsavg.id.avg)
tail(gddsaSFNormcountsavg.id.avg)
gddsaSFNormcountsavg.selected_genes = merge(gddsaSFNormcountsavg.id.avg, selectedgenes.cords, by="id", all.x=FALSE)
head(gddsaSFNormcountsavg.selected_genes)
gddsaSFNormcountsavg.selected_genes_sort <- gddsaSFNormcountsavg.selected_genes[order(gddsaSFNormcountsavg.selected_genes$chr, gddsaSFNormcountsavg.selected_genes$start),]
head(gddsaSFNormcountsavg.selected_genes_sort)
dim(gddsaSFNormcountsavg.selected_genes_sort)
gddsaSFNormcountsavg.selected_genes_sort <- gddsaSFNormcountsavg.selected_genes_sort[c(1,13:34, 2:12),]
head(gddsaSFNormcountsavg.selected_genes_sort)
gddsaSFNormcountsavg.selected_genes_sort_heat1 <- gddsaSFNormcountsavg.selected_genes_sort[,c(9,2:5)]
rownames(gddsaSFNormcountsavg.selected_genes_sort_heat1)
rownames(gddsaSFNormcountsavg.selected_genes_sort_heat1)=gddsaSFNormcountsavg.selected_genes_sort_heat1[,1]
rownames(gddsaSFNormcountsavg.selected_genes_sort_heat1)
head(gddsaSFNormcountsavg.selected_genes_sort_heat1)
write.table(gddsaSFNormcountsavg.selected_genes_sort_heat1, "gddsaSFNormcountsavg.selected_genes_sort_heat1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb")

gddsaSFNormcountsavg.selected_genes_sort_heat1 <- read.table("gddsaSFNormcountsavg.selected_genes_sort_heat1.txt", header = T)
head(gddsaSFNormcountsavg.selected_genes_sort_heat1)
gddsaSFNormcountsavg.selected_genes_sort_heat1["WTratio"] <- gddsaSFNormcountsavg.selected_genes_sort_heat1$JB1_WTJF1 / (gddsaSFNormcountsavg.selected_genes_sort_heat1$JB1_WTB6 + gddsaSFNormcountsavg.selected_genes_sort_heat1$JB1_WTJF1)
gddsaSFNormcountsavg.selected_genes_sort_heat1["ZFP57KOratio"]  <- gddsaSFNormcountsavg.selected_genes_sort_heat1$JB1_ZFP57KOJF1 / (gddsaSFNormcountsavg.selected_genes_sort_heat1$JB1_ZFP57KOB6 + gddsaSFNormcountsavg.selected_genes_sort_heat1$JB1_ZFP57KOJF1)
gddsaSFNormcountsavg.selected_genes_sort_heat1["B6ratio"] <- log2(gddsaSFNormcountsavg.selected_genes_sort_heat1$JB1_ZFP57KOB6 / gddsaSFNormcountsavg.selected_genes_sort_heat1$JB1_WTB6)
gddsaSFNormcountsavg.selected_genes_sort_heat1["JF1ratio"]  <- log2(gddsaSFNormcountsavg.selected_genes_sort_heat1$JB1_ZFP57KOJF1 / gddsaSFNormcountsavg.selected_genes_sort_heat1$JB1_WTJF1)
write.table(gddsaSFNormcountsavg.selected_genes_sort_heat1, "gddsaSFNormcountsavg.selected_genes_sort_heat1.withratio.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)

head(gddsaSFNormcountsavg.selected_genes_sort_heat1)
library(pheatmap)
library(RColorBrewer)
#Predict Mnoallecity or Biallecity
#pheatmap(data1,treeheight_row = 0, cluster_cols=F, cluster_rows=F, treeheight_col = 0, gaps_col =NULL, gaps_row = NULL, border_color = "black", breaks = breaksList,color= colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)))
breaksList1 = c(0, 0.31, 0.67)
pheatmap(gddsaSFNormcountsavg.selected_genes_sort_heat1[,6:7],
         color = colorRampPalette(c("#00B0F0","#f8de38","#f838b0"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "black",
         cluster_cols=F, 
         cluster_rows=F,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12,
         filename = "z_gddsaSFNormcountsavg.selected_genes_sort_heat1_Interpretation.png")
dev.off()
#Predict silent or active allele
 
#Copy gddsaSFNormcountsavg.selected_genes_sort_heat1.withratio.txt file in excel and separate silent and active allele and save as gddsaSFNormcountsavg.selected_genes_sort_heat1.withratio.reorder.txt. The excel sheet is in 5mb folder names as order heatmap silent active.xlsx
gddsaSFNormcountsavg.selected_genes_sort_heat1.withratio.reorder <- read.table("gddsaSFNormcountsavg.selected_genes_sort_heat1.withratio.reorder.txt", header = T)
head(gddsaSFNormcountsavg.selected_genes_sort_heat1.withratio.reorder)
rownames(gddsaSFNormcountsavg.selected_genes_sort_heat1.withratio.reorder) <- gddsaSFNormcountsavg.selected_genes_sort_heat1.withratio.reorder[,1]
gddsaSFNormcountsavg.selected_genes_sort_heat1.withratio.reorder <- gddsaSFNormcountsavg.selected_genes_sort_heat1.withratio.reorder[,-1]
gddsaSFNormcountsavg.selected_genes_sort_heat1.withratio.reorder
breaksList1 = c(-300,-200,-0.15,0.15)
pheatmap(gddsaSFNormcountsavg.selected_genes_sort_heat1.withratio.reorder[,1:2],
         color = colorRampPalette(c("darkgrey","#a10517","white","#35b04b"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "black",
         cluster_cols=F, 
         cluster_rows=F,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12,
         filename = "z_gddsaSFNormcountsavg.selected_genes_sort_heat1_Interpretation_updown.png")

dev.off()
#Heatmap for imprinted genes
grep ENS imprinted_gene_rev.sort.txt > imprinted_gene_cords.bed #grep coordinates

imprinted_gene_cords <- read.table("imprinted_gene_cords.bed", header = F)
head(imprinted_gene_cords)
colnames(imprinted_gene_cords) <- c("chr", "start", "end", "id", "Gene")
head(imprinted_gene_cords)
dim(imprinted_gene_cords)
head(z_gddsaSFNormcountsavg)
dim(z_gddsaSFNormcountsavg)
z_gddsaSFNormcountsavg_sort <- data.frame(cbind(rownames(z_gddsaSFNormcountsavg), z_gddsaSFNormcountsavg))
head(z_gddsaSFNormcountsavg_sort)
colnames(z_gddsaSFNormcountsavg_sort) <- c("id", "JB1_WTB6", "JB1_WTJF1", "JB1_ZFP57KOB6", "JB1_ZFP57KOJF1")
z_gddsaSFNormcountsavg_imprinted_sort = merge(z_gddsaSFNormcountsavg_sort, imprinted_gene_cords, by="id", all.x=F)
head(z_gddsaSFNormcountsavg_imprinted_sort)
dim(z_gddsaSFNormcountsavg_imprinted_sort)
write.table(z_gddsaSFNormcountsavg_imprinted_sort, "z_gddsaSFNormcountsavg_imprinted_sort.txt", quote = F, append = F )
library(ggplot2)
library(gplots)
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb")
z_gddsaSFNormcountsavg_imprinted_sortre <- z_gddsaSFNormcountsavg_imprinted_sort[,c(6:9,1:5)]
head(z_gddsaSFNormcountsavg_imprinted_sortre)
write.table(z_gddsaSFNormcountsavg_imprinted_sortre, "z_gddsaSFNormcountsavg_imprinted_sortre.txt", quote = F, append= F, row.names = F)
z_gddsa_avgimppos <- read.table("z_gddsaSFNormcountsavg_imprinted_sortre.txt",header=T)
head(z_gddsa_avgimppos)
rownames(z_gddsa_avgimppos) <- z_gddsa_avgimppos$Gene
z_gddsa_avgimppos1 = z_gddsa_avgimppos[,6:9]
z_gddsa_avgimppos2 = as.matrix(z_gddsa_avgimppos1)
head(z_gddsa_avgimppos2)
dim(z_gddsa_avgimppos2)
colfunc <- colorRampPalette(c("navy","white", "#C71111"))
#heatmap.2(z_gddsa_avgimp1, Colv = "NA",trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=1, font=3, cexCol = 1, margins =c(3,8), srtCol = 45, breaks = seq(-5,5, length.out = 100))
heatmap.2(z_gddsa_avgimppos2, Colv = "NA",Rowv ='NA',trace = "none", col = colfunc ,
          lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5),
          density.info=c("none"), dendrogram="none", scale = "none", 
          sepwidth=c(0.001, 0.001), cexRow=1, font=3, cexCol = 0.8, 
          sepcolor="black", margins =c(6,8), srtCol = 45, 
          breaks = seq(-2,2, length.out = 100), colsep=1:ncol(z_gddsa_avgimppos2),
          rowsep=1:nrow(z_gddsa_avgimppos2))

#save as z_gddsa_avgimppos2.svg
# load package
library(pheatmap)
library(RColorBrewer)
breaksList = seq(-3, 3)

pheatmap(z_gddsa_avgimppos2,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)),
         breaks = breaksList,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 1,
         cellwidth = 12, 
         cellheight = 12,
         filename = "z_gddsa_avgimppos2.png")

#save as z_gddsa_avgimppos2.svg

#Fisher test and tag chromsome positions
gddsaSFNormcountsavg_chr_gene1 <- read.table("gddsaSFNormcountsavg_chr_gene1.txt", header = T)
head(gddsaSFNormcountsavg_chr_gene1)
dim(gddsaSFNormcountsavg_chr_gene1)
rownames(gddsaSFNormcountsavg_chr_gene1) <- paste(gddsaSFNormcountsavg_chr_gene1$Genes,gddsaSFNormcountsavg_chr_gene1$id, sep = "%")
gddsaSFNormcountsavg_chr_gene1_fisher <- gddsaSFNormcountsavg_chr_gene1[,6:9]
head(gddsaSFNormcountsavg_chr_gene1_fisher)
gddsaSFNormcountsavg_chr_gene1_fisher1 <- apply(gddsaSFNormcountsavg_chr_gene1_fisher,1, function(x) fisher.test(matrix(x,nrow =2))$p.value)
#warning will be observed because we need only integers for runnung this function.
#In fisher.test(matrix(x, nrow = 2)) :'x' has been rounded to integer: Mean relative difference:
head(gddsaSFNormcountsavg_chr_gene1_fisher1)
gddsaSFNormcountsavg_chr_gene1_fisher["fisherpvalue"] <- data.frame(gddsaSFNormcountsavg_chr_gene1_fisher1)
head(gddsaSFNormcountsavg_chr_gene1_fisher)
gddsaSFNormcountsavg_chr_gene1_fisher["fisheradjpval"] <- p.adjust(gddsaSFNormcountsavg_chr_gene1_fisher1,method="BH")
head(gddsaSFNormcountsavg_chr_gene1_fisher)
write.table(gddsaSFNormcountsavg_chr_gene1_fisher, "gddsaSFNormcountsavg_chr_gene1_fisher.txt", sep="\t", quote = FALSE, append = FALSE)

gddsaSFNormcountsavg_chr_gene1_fishercoordinate <- cbind.data.frame(gddsaSFNormcountsavg_chr_gene1, gddsaSFNormcountsavg_chr_gene1_fisher)
head(gddsaSFNormcountsavg_chr_gene1_fishercoordinate)
write.table(gddsaSFNormcountsavg_chr_gene1_fishercoordinate, "gddsaSFNormcountsavg_chr_gene1_fishercoordinate.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)
writexl::write_xlsx(gddsaSFNormcountsavg_chr_gene1_fishercoordinate, "JB1_NPCs_T12_gddsaSFNormcountsavg_chr_gene1_fishercoordinate_mm10.xlsx")

gddsaSFNormcounts_chr_gene1_fishercoordinate <- cbind.data.frame(gddsaSFNormcounts_chr_gene1, gddsaSFNormcountsavg_chr_gene1_fisher)
head(gddsaSFNormcounts_chr_gene1_fishercoordinate)
tail(gddsaSFNormcounts_chr_gene1_fishercoordinate)
write.table(gddsaSFNormcounts_chr_gene1_fishercoordinate, "gddsaSFNormcounts_chr_gene1_fishercoordinate.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)
writexl::write_xlsx(gddsaSFNormcounts_chr_gene1_fishercoordinate, "JB1_NPCs_T12_gddsaSFNormcounts_chr_gene1_fishercoordinate_mm10.xlsx")

#WithID
gddsaSFNormcountsavg_chr_gene1_fisher["GeneID"] <- rownames(gddsaSFNormcountsavg_chr_gene1_fisher)
gddsaSFNormcountsavg_chr_gene1_fisher_0.05 <- gddsaSFNormcountsavg_chr_gene1_fisher[which(gddsaSFNormcountsavg_chr_gene1_fisher$fisheradjpval <0.05),]
head(gddsaSFNormcountsavg_chr_gene1_fisher_0.05)
dim(gddsaSFNormcountsavg_chr_gene1_fisher_0.05)
write.table(gddsaSFNormcountsavg_chr_gene1_fisher_0.05, "gddsaSFNormcountsavg_chr_gene1_fisher_0.05.txt", sep="\t", quote = FALSE, append = FALSE)

#WithID
gddsaSFNormcountsavg_chr_gene1_fisher_0.1 <- gddsaSFNormcountsavg_chr_gene1_fisher[which(gddsaSFNormcountsavg_chr_gene1_fisher$fisheradjpval <0.1),]
head(gddsaSFNormcountsavg_chr_gene1_fisher_0.1)
dim(gddsaSFNormcountsavg_chr_gene1_fisher_0.1)
write.table(gddsaSFNormcountsavg_chr_gene1_fisher_0.1, "gddsaSFNormcountsavg_chr_gene1_fisher_0.1.txt", sep="\t", quote = FALSE, append = FALSE)


#All imprinted genes maternal ratio
head(gddsaSFNormcountsavg_chr_gene1)
gddsaSFNormcountsavg_chr_gene1_imprinted_sort = merge(gddsaSFNormcountsavg_chr_gene1, imprinted_gene_cords, by="id", all.x=F)
head(gddsaSFNormcountsavg_chr_gene1_imprinted_sort)
dim(gddsaSFNormcountsavg_chr_gene1_imprinted_sort)
gddsaSFNormcountsavg_chr_gene1_imprinted_sortre <- gddsaSFNormcountsavg_chr_gene1_imprinted_sort[,c(2:4,1,5:9)]
head(gddsaSFNormcountsavg_chr_gene1_imprinted_sortre)
gddsaSFNormcountsavg_chr_gene1_imprinted_sortre["WTmatratio"] <- gddsaSFNormcountsavg_chr_gene1_imprinted_sortre$JB1_WTJF1 / (gddsaSFNormcountsavg_chr_gene1_imprinted_sortre$JB1_WTB6 + gddsaSFNormcountsavg_chr_gene1_imprinted_sortre$JB1_WTJF1)
gddsaSFNormcountsavg_chr_gene1_imprinted_sortre["Zfp57komatratio"] <- gddsaSFNormcountsavg_chr_gene1_imprinted_sortre$JB1_ZFP57KOJF1 / (gddsaSFNormcountsavg_chr_gene1_imprinted_sortre$JB1_ZFP57KOB6 + gddsaSFNormcountsavg_chr_gene1_imprinted_sortre$JB1_ZFP57KOJF1)
head(gddsaSFNormcountsavg_chr_gene1_imprinted_sortre)
write.table(gddsaSFNormcountsavg_chr_gene1_imprinted_sortre, "gddsaSFNormcountsavg_chr_gene1_imprinted_sortre.txt", quote = F, append = F, row.names = F)

#grep chrX gddsaSFNormcountsavg_chr_gene1_imprinted_sortre.txt -v | awk '{if($10 <= 0.31 || $10 >= 0.67) print $0}' > gddsaSFNormcountsavg_chr_gene1_imprinted_sortre.filt.txt

#Scatter Plot
head(gddsaSFNormcountsavg_chr_gene1)
gddsaSFNormcountsavg_chr_gene1_ratio <- gddsaSFNormcountsavg_chr_gene1
gddsaSFNormcountsavg_chr_gene1_ratio["matWTratio"] <- gddsaSFNormcountsavg_chr_gene1_ratio$JB1_WTJF1/(gddsaSFNormcountsavg_chr_gene1_ratio$JB1_WTB6 + gddsaSFNormcountsavg_chr_gene1_ratio$JB1_WTJF1)
head(gddsaSFNormcountsavg_chr_gene1_ratio)
gddsaSFNormcountsavg_chr_gene1_ratio["matZFP57KOratio"] <- gddsaSFNormcountsavg_chr_gene1_ratio$JB1_ZFP57KOJF1/(gddsaSFNormcountsavg_chr_gene1_ratio$JB1_ZFP57KOB6 + gddsaSFNormcountsavg_chr_gene1_ratio$JB1_ZFP57KOJF1)
head(gddsaSFNormcountsavg_chr_gene1_ratio)

plot(gddsaSFNormcountsavg_chr_gene1_ratio$matWTratio, col="darkgreen", pch =16, cex = 0.3)
par(new=T)
plot(gddsaSFNormcountsavg_chr_gene1_ratio$matZFP57KOratio, col="darkred", pch =16, cex = 0.3)

gddsaSFNormcountsavg_chr_gene1_ratioR <- data.frame(gddsaSFNormcountsavg_chr_gene1_ratio[,10:11])
head(gddsaSFNormcountsavg_chr_gene1_ratioR[c(4500:4900),])
gddsaSFNormcountsavgratioRange <- gddsaSFNormcountsavg_chr_gene1_ratioR[c(4500:4900),]
head(gddsaSFNormcountsavgratioRange)
plot(gddsaSFNormcountsavgratioRange$WTratio, col="darkgreen", pch =16, cex = 0.8)
par(new=T)
plot(gddsaSFNormcountsavgratioRange$ZFP57KOratio, col="darkred", pch =16, cex = 0.8)


#Basilia provided the an excel sheet file which I saved as Table_for_ankit _27_july.csv. 
#I created a same sheet: fgrep -f basilia_genes_for_dotplot.txt gddsaSFNormcounts_chr_gene1.txt -w | sort -k4,4 > basilia_genes_for_dotplot_gddsaSFNormcounts.txt

#Manually done
#Open basilia_genes_for_dotplot_gddsaSFNormcounts.txt and if 0 occurs in both allele, we replaced it with 1 (atleast one read), the other allele was still high
#save as basilia_genes_for_dotplot_gddsaSFNormcounts_re.txt
basilia_genes_for_dotplot_gddsaSFNormcounts <- read.table("basilia_genes_for_dotplot_gddsaSFNormcounts_re.txt", header = F)
colnames(basilia_genes_for_dotplot_gddsaSFNormcounts) <- colnames(gddsaSFNormcounts_chr_gene1)
head(basilia_genes_for_dotplot_gddsaSFNormcounts)
basilia_genes_for_dotplot_gddsaSFNormcounts["WT_allele_ratio_1"] <- basilia_genes_for_dotplot_gddsaSFNormcounts$JB1_WT_Rep1.JF1 / (basilia_genes_for_dotplot_gddsaSFNormcounts$JB1_WT_Rep1.B6 + basilia_genes_for_dotplot_gddsaSFNormcounts$JB1_WT_Rep1.JF1)
basilia_genes_for_dotplot_gddsaSFNormcounts["WT_allele_ratio_2"] <- basilia_genes_for_dotplot_gddsaSFNormcounts$JB1_WT_Rep2.JF1 / (basilia_genes_for_dotplot_gddsaSFNormcounts$JB1_WT_Rep2.B6 + basilia_genes_for_dotplot_gddsaSFNormcounts$JB1_WT_Rep2.JF1)
basilia_genes_for_dotplot_gddsaSFNormcounts["WT_allele_ratio"] <- (basilia_genes_for_dotplot_gddsaSFNormcounts$WT_allele_ratio_1 + basilia_genes_for_dotplot_gddsaSFNormcounts$WT_allele_ratio_2)/2
basilia_genes_for_dotplot_gddsaSFNormcounts["ZFP57KO_allele_ratio_1"] <- basilia_genes_for_dotplot_gddsaSFNormcounts$JB1_ZFP57_KO_Rep1.JF1 / (basilia_genes_for_dotplot_gddsaSFNormcounts$JB1_ZFP57_KO_Rep1.B6 + basilia_genes_for_dotplot_gddsaSFNormcounts$JB1_ZFP57_KO_Rep1.JF1)
basilia_genes_for_dotplot_gddsaSFNormcounts["ZFP57KO_allele_ratio_2"] <- basilia_genes_for_dotplot_gddsaSFNormcounts$JB1_ZFP57_KO_Rep2.JF1 / (basilia_genes_for_dotplot_gddsaSFNormcounts$JB1_ZFP57_KO_Rep2.B6 + basilia_genes_for_dotplot_gddsaSFNormcounts$JB1_ZFP57_KO_Rep2.JF1)
basilia_genes_for_dotplot_gddsaSFNormcounts["ZFP57KO_allele_ratio"] <- (basilia_genes_for_dotplot_gddsaSFNormcounts$ZFP57KO_allele_ratio_1 + basilia_genes_for_dotplot_gddsaSFNormcounts$ZFP57KO_allele_ratio_2)/2

data_imp_gene_indiv_asp_norm <- basilia_genes_for_dotplot_gddsaSFNormcounts
# order
data_imp_gene_indiv_asp_norm <- data_imp_gene_indiv_asp_norm[order(data_imp_gene_indiv_asp_norm$WT_allele_ratio),]
head(data_imp_gene_indiv_asp_norm)
dim(data_imp_gene_indiv_asp_norm)
write.table(data_imp_gene_indiv_asp_norm, "data_imp_gene_indiv_asp_norm.txt", sep = "\t", quote = F, append = F, row.names = F)

summary(data_imp_gene_indiv_asp_norm[,c(16,19)])
data_imp_gene_indiv_asp_norm_re <- data_imp_gene_indiv_asp_norm[,c(4,16,19)]
head(data_imp_gene_indiv_asp_norm_re)
colnames(data_imp_gene_indiv_asp_norm_re) <- c("Gene", "WT_avg_allele_ratio", "ZFP57KO_avg_allele_ratio")
rownames(data_imp_gene_indiv_asp_norm_re) <- data_imp_gene_indiv_asp_norm_re[,1]
data_imp_gene_indiv_asp_norm_rebarplot <- data_imp_gene_indiv_asp_norm_re[,-1]
data_imp_gene_indiv_asp_norm_rebarplot <- as.matrix(t(data_imp_gene_indiv_asp_norm_rebarplot))
head(data_imp_gene_indiv_asp_norm_rebarplot)
dim(data_imp_gene_indiv_asp_norm_rebarplot)
#barplot(data_imp_gene_indiv_asp_norm_rebarplot, beside=TRUE, horiz=TRUE, col=c("red","blue","orange", "green"), xpd=FALSE,  xlim = c(0,1))
data_imp_gene_indiv_asp_norm_rebarplot1 <- stack(data_imp_gene_indiv_asp_norm_rebarplot)
head(data_imp_gene_indiv_asp_norm_rebarplot1)

colnames(data_imp_gene_indiv_asp_norm_rebarplot1) <- c("Group", "Gene", "Ratio")
data_imp_gene_indiv_asp_norm_rebarplot1 <- data.frame(data_imp_gene_indiv_asp_norm_rebarplot1)
head(data_imp_gene_indiv_asp_norm_rebarplot1)
str(data_imp_gene_indiv_asp_norm_rebarplot1)
data_imp_gene_indiv_asp_norm_rebarplot1_half <- data_imp_gene_indiv_asp_norm_rebarplot1
data_imp_gene_indiv_asp_norm_rebarplot1_half["Ratio"] <- as.numeric(data_imp_gene_indiv_asp_norm_rebarplot1_half$Ratio) - 0.5
head(data_imp_gene_indiv_asp_norm_rebarplot1_half)
write.table(data_imp_gene_indiv_asp_norm_rebarplot1_half, "data_imp_gene_indiv_asp_norm_rebarplot1_half.txt", sep = "\t", quote = F, append = F)
library(ggplot2)
# Basic barplot

#way1 : ggplot2
pbar <-ggplot(data=data_imp_gene_indiv_asp_norm_rebarplot1_half, aes(x=Gene, y=Ratio, label=Group, color=Group, fill=Group)) +
  geom_bar(width = 0.5,stat="identity",  position=position_dodge(), size=0.6)
pbar + coord_flip()+
  scale_fill_manual(values=c("darkgreen","darkred"))+
  scale_color_manual(values=c("darkgreen","darkred")) + theme_bw()
ggsave("barplot_data_imp_gene_indiv_asp_norm_rebarplot1_half.svg", width=13, height=25, units="cm", dpi=96)
ggsave("barplot_data_imp_gene_indiv_asp_norm_rebarplot1_half.jpg", width=13, height=25, units="cm", dpi=96)


#way2 : ggpubr
library(ggpubr)
ggbarplot(data_imp_gene_indiv_asp_norm_rebarplot1_half, x = "Gene", y = "Ratio", fill = "Group", color = "white", palette = c("darkgreen","darkred"),sort.by.groups = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Ratio",
          legend.title = "Expression",
          rotate = TRUE, position = position_dodge(),
          ggtheme = theme_classic())
ggbarplot(data_imp_gene_indiv_asp_norm_rebarplot1_half, x = "Gene", y = "Ratio", fill = "Group", color = "white", palette = c("darkgreen","darkred"),sort.by.groups = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Ratio",
          legend.title = "Expression",
          rotate = TRUE, position = position_dodge(),
          ggtheme = theme_classic())

head(data_imp_gene_indiv_asp_norm)
data_imp_gene_indiv_asp_norm_are <- data_imp_gene_indiv_asp_norm[,c(4,14:15,17:18)]
head(data_imp_gene_indiv_asp_norm_are)
data_imp_gene_indiv_asp_norm_sortre <- data_imp_gene_indiv_asp_norm_are
head(data_imp_gene_indiv_asp_norm_sortre)
rownames(data_imp_gene_indiv_asp_norm_sortre) <- data_imp_gene_indiv_asp_norm_sortre[,1]
data_imp_gene_indiv_asp_norm_sortrebarplot <- data_imp_gene_indiv_asp_norm_sortre[,-1]
data_imp_gene_indiv_asp_norm_sortrebarplot <- as.matrix(t(data_imp_gene_indiv_asp_norm_sortrebarplot))
head(data_imp_gene_indiv_asp_norm_sortrebarplot)
dim(data_imp_gene_indiv_asp_norm_sortrebarplot)
#barplot(data_imp_gene_indiv_asp_norm_sortrebarplot, beside=TRUE, horiz=TRUE, col=c("red","blue","orange", "green"), xpd=FALSE,  xlim = c(0,1))


data_imp_gene_indiv_asp_norm_sortrebarplot1 <- stack(data_imp_gene_indiv_asp_norm_sortrebarplot)
head(data_imp_gene_indiv_asp_norm_sortrebarplot1)
colnames(data_imp_gene_indiv_asp_norm_sortrebarplot1) <- c("Group", "Gene", "Ratio")
data_imp_gene_indiv_asp_norm_sortrebarplot1 <- data.frame(data_imp_gene_indiv_asp_norm_sortrebarplot1)
head(data_imp_gene_indiv_asp_norm_sortrebarplot1)
str(data_imp_gene_indiv_asp_norm_sortrebarplot1)
data_imp_gene_indiv_asp_norm_sortrebarplot1_half <- data_imp_gene_indiv_asp_norm_sortrebarplot1
data_imp_gene_indiv_asp_norm_sortrebarplot1_half["Ratio"] <- as.numeric(data_imp_gene_indiv_asp_norm_sortrebarplot1_half$Ratio) - 0.5
head(data_imp_gene_indiv_asp_norm_sortrebarplot1_half)
library(ggplot2)
# Basic barplot

#way1 : ggplot2
pbar <-ggplot(data=data_imp_gene_indiv_asp_norm_sortrebarplot1_half, aes(x=Gene, y=Ratio, label=Group, color=Group, fill=Group)) +
  geom_bar(width = 0.5,stat="identity",  position=position_dodge(), size=0.6)
pbar + coord_flip()+
  scale_fill_manual(values=c("darkgreen","darkgreen","darkred","darkred"))+
  scale_color_manual(values=c("darkgreen","darkgreen","darkred","darkred")) + theme_bw()
ggsave("data_imp_gene_indiv_asp_norm_sortrebarplot1_half.svg", width=13, height=25, units="cm", dpi=96)
ggsave("data_imp_gene_indiv_asp_norm_sortrebarplot1_half.jpg", width=13, height=25, units="cm", dpi=96)
#theme(legend.position = "bottom", legend.direction = "horizontal", panel.background = element_blank(), axis.line.x = element_line(size = .6, colour = "black"))


pdot <-ggplot(data=data_imp_gene_indiv_asp_norm_sortrebarplot1_half, aes(x=Gene, y=Ratio, label=Group, color="black", fill=Group)) + 
  geom_hline(yintercept = c(-0.17,0,0.17), colour = "grey", linetype=c("dashed","solid","dashed")) + 
  geom_dotplot(binaxis='y', stackdir='centerwhole', stackratio=1.5, dotsize=1.2) + theme_classic() 
pdot + coord_flip()+ ylim(c(-0.5,0.5))+
  scale_fill_manual(values=c("green","green","red","red"))+
  scale_color_manual(values=c("black","black","black","black")) 

ggsave("dotplot_data_imp_gene_indiv_asp_norm_sortrebarplot1_half.svg", width=13, height=24, units="cm", dpi=96)
ggsave("dotplot_data_imp_gene_indiv_asp_norm_sortrebarplot1_half.jpg", width=13, height=24, units="cm", dpi=96)

pdotsquare <-ggplot(data=data_imp_gene_indiv_asp_norm_sortrebarplot1_half, aes(x=Gene, y=Ratio, label=Group, color="black", fill=Group)) + 
  geom_hline(yintercept = c(-0.17,0,0.17), colour = "grey", linetype=c("dashed","solid","dashed")) + 
  geom_point(aes(shape=Group, color=Group, fill=Group, size=Group), position = position_dodge(width = 0.50)) + theme_classic()  

pdotsquare + coord_flip()+ ylim(c(-0.5,0.5))+
  scale_fill_manual(values=c("green","green","red","red"))+
  scale_color_manual(values=c("black","black","black","black")) +
  scale_shape_manual(values=c(22,24,22,24))+
  scale_size_manual(values=c(1.5,1.5,1.5,1.5)) 

ggsave("trisquareplot_data_imp_gene_indiv_asp_norm_sortrebarplot1_halfdodge.svg", width=13, height=24, units="cm", dpi=96)
ggsave("trisquareplot_data_imp_gene_indiv_asp_norm_sortrebarplot1_halfdodge.jpg", width=13, height=24, units="cm", dpi=96)

pdotsquare <-ggplot(data=data_imp_gene_indiv_asp_norm_sortrebarplot1_half, aes(x=Gene, y=Ratio, label=Group, color="black", fill=Group)) + 
  geom_hline(yintercept = c(-0.17,0,0.17), colour = "grey", linetype=c("dashed","solid","dashed")) + 
  geom_point(aes(shape=Group, color=Group, fill=Group, size=Group), position = position_dodge(width = 0.1)) + theme_classic()  

pdotsquare + coord_flip()+ ylim(c(-0.5,0.5))+
  scale_fill_manual(values=c("green","green","red","red"))+
  scale_color_manual(values=c("black","black","black","black")) +
  scale_shape_manual(values=c(22,24,22,24))+
  scale_size_manual(values=c(1.5,1.5,1.5,1.5)) 

ggsave("trisquareplot_data_imp_gene_indiv_asp_norm_sortrebarplot1_half.svg", width=13, height=24, units="cm", dpi=96)
ggsave("trisquareplot_data_imp_gene_indiv_asp_norm_sortrebarplot1_half.jpg", width=13, height=24, units="cm", dpi=96)


#way2 : ggpubr
#library(ggpubr)
#ggbarplot(data_imp_gene_indiv_asp_norm_sortrebarplot1_half, x = "Gene", y = "Ratio", fill = "Group", color = "white", palette = c("darkgreen","darkgreen","darkred","darkred"),sort.by.groups = FALSE,
#          x.text.angle = 90,          # Rotate vertically x axis texts
#          ylab = "Ratio",
#          legend.title = "Expression",
#          rotate = TRUE, position = position_dodge(),
#          ggtheme = theme_classic())
#ggdotplot(data_imp_gene_indiv_asp_norm_sortrebarplot1_half, x = "Gene", y = "Ratio", fill = "Group", color = " black", palette = c("darkgreen","darkgreen","darkred","darkred"),sort.by.groups = FALSE,
#          x.text.angle = 90,          # Rotate vertically x axis texts
#          ylab = "Ratio",
#          merge = TRUE,
#          size=0.5,
#          binwidth = 0.05,
#          legend.title = "Expression",
#          rotate = TRUE,
#          repel = 2,ggsave("dotplot_data_imp_gene_indiv_asp_norm_sortrebarplot1_half.jpg", width=13, height=30, units="cm", dpi=96)
#          ggtheme =  theme_classic()) + theme_cleveland() + geom_hline(yintercept = 0.0, colour = "grey") 
#ggsave("dotplot_themeclassic_data_imp_gene_indiv_asp_norm_sortrebarplot1_half.svg", width=13, height=25, units="cm", dpi=96)
#theme(legend.position = "bottom", legend.direction = "horizontal", panel.background = element_blank(), axis.line.x = element_line(size = .6, colour = "black"))
#fgrep -f selected_imprinted_genes.txt /home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20_rearranged.txt -w | grep Mrpl23-ps1 -v > selected_imprinted_genes.cords.txt

#Cluster Dotplot of 39 imprinted genes. 
#Basilia set of 37 genes  and Blcap , Adam23 added laterand more later Trappc9 on so 40
basilia_genes_for_clusterdotplot_as_bulk <- read.table("basilia_genes_for_clusterdotplot_as_bulk.txt", header = F)
colnames(basilia_genes_for_clusterdotplot_as_bulk) <- c("Genes")
head(gddsaSFNormcounts_chr_gene1)
dim(gddsaSFNormcounts_chr_gene1)
gddsaSFNormcounts_chr_gene1_clusterdotplotgenes <- merge(basilia_genes_for_clusterdotplot_as_bulk,gddsaSFNormcounts_chr_gene1, by ="Genes", all.x=FALSE,sort = F)
head(gddsaSFNormcounts_chr_gene1_clusterdotplotgenes)
dim(gddsaSFNormcounts_chr_gene1_clusterdotplotgenes)

write.csv(gddsaSFNormcounts_chr_gene1_clusterdotplotgenes, "gddsaSFNormcounts_chr_gene1_clusterdotplotgenes.csv", row.names = F)


#Manually done
#Open gddsaSFNormcounts_chr_gene1_clusterdotplotgenes.csv and if 0 occurs in both allele, we replaced it with 1 (atleast one read), the other allele was still high
#save as gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re.txt, do not copy Genes word so that header can kept true


gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re <- read.table("gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re.txt", header = T)
gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re["WT_allele_ratio_1"] <- gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re$JB1_WT_Rep1.JF1 / (gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re$JB1_WT_Rep1.B6 + gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re$JB1_WT_Rep1.JF1)
gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re["WT_allele_ratio_2"] <- gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re$JB1_WT_Rep2.JF1 / (gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re$JB1_WT_Rep2.B6 + gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re$JB1_WT_Rep2.JF1)
gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re["WT_allele_ratio"] <- (gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re$WT_allele_ratio_1 + gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re$WT_allele_ratio_2)/2
gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re["ZFP57KO_allele_ratio_1"] <- gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re$JB1_ZFP57_KO_Rep1.JF1 / (gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re$JB1_ZFP57_KO_Rep1.B6 + gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re$JB1_ZFP57_KO_Rep1.JF1)
gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re["ZFP57KO_allele_ratio_2"] <- gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re$JB1_ZFP57_KO_Rep2.JF1 / (gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re$JB1_ZFP57_KO_Rep2.B6 + gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re$JB1_ZFP57_KO_Rep2.JF1)
gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re["ZFP57KO_allele_ratio"] <- (gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re$ZFP57KO_allele_ratio_1 + gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re$ZFP57KO_allele_ratio_2)/2

head(gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re,2)
cluster_imp_gene_indiv_asp_norm <- gddsaSFNormcounts_chr_gene1_clusterdotplotgenes_re

head(cluster_imp_gene_indiv_asp_norm)
dim(cluster_imp_gene_indiv_asp_norm)
write.table(cluster_imp_gene_indiv_asp_norm, "cluster_imp_gene_indiv_asp_norm.txt", sep = "\t", quote = F, append = F, row.names = T)

summary(cluster_imp_gene_indiv_asp_norm[,c(13:14,16:17)])

#Separate replicates
cluster_imp_gene_indiv_asp_norm_are <- cluster_imp_gene_indiv_asp_norm[,c(13:14,16:17)]
head(cluster_imp_gene_indiv_asp_norm_are)
cluster_imp_gene_indiv_asp_norm_sortre <- cluster_imp_gene_indiv_asp_norm_are
head(cluster_imp_gene_indiv_asp_norm_sortre)
cluster_imp_gene_indiv_asp_norm_sortrebarplot <- cluster_imp_gene_indiv_asp_norm_sortre
cluster_imp_gene_indiv_asp_norm_sortrebarplot <- as.matrix(t(cluster_imp_gene_indiv_asp_norm_sortrebarplot))
head(cluster_imp_gene_indiv_asp_norm_sortrebarplot)
dim(cluster_imp_gene_indiv_asp_norm_sortrebarplot)
cluster_imp_gene_indiv_asp_norm_sortrebarplot1 <- stack(cluster_imp_gene_indiv_asp_norm_sortrebarplot)
head(cluster_imp_gene_indiv_asp_norm_sortrebarplot1)
colnames(cluster_imp_gene_indiv_asp_norm_sortrebarplot1) <- c("Group", "Gene", "Ratio")
cluster_imp_gene_indiv_asp_norm_sortrebarplot1 <- data.frame(cluster_imp_gene_indiv_asp_norm_sortrebarplot1)
head(cluster_imp_gene_indiv_asp_norm_sortrebarplot1)
dim(cluster_imp_gene_indiv_asp_norm_sortrebarplot1)
cluster_imp_gene_indiv_asp_norm_sortrebarplot1["Group.1"] <- data.frame(rep(c(rep("WT_allele_ratio",2),rep("ZFP57KO_allele_ratio",2)), 40))

head(cluster_imp_gene_indiv_asp_norm_sortrebarplot1)
dim(cluster_imp_gene_indiv_asp_norm_sortrebarplot1)
str(cluster_imp_gene_indiv_asp_norm_sortrebarplot1)
cluster_imp_gene_indiv_asp_norm_sortrebarplot1_half <- cluster_imp_gene_indiv_asp_norm_sortrebarplot1
cluster_imp_gene_indiv_asp_norm_sortrebarplot1_half["Ratio"] <- as.numeric(cluster_imp_gene_indiv_asp_norm_sortrebarplot1_half$Ratio) - 0.5
head(cluster_imp_gene_indiv_asp_norm_sortrebarplot1_half)
library(ggplot2)
# Dotplot for Cluster
# Change the position
pcluster <-ggplot(cluster_imp_gene_indiv_asp_norm_sortrebarplot1_half, aes(x=Gene, y=Ratio, fill= Group.1)) + geom_hline(yintercept = c(-0.17,0,0.17), colour = "grey", linetype=c("dashed","solid","dashed"))+
  geom_boxplot(position=position_dodge(0.8))+
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.8))
pcluster + scale_fill_manual(values=c("green","red"))+ylim(c(-0.5,0.5))+
  scale_color_manual(values=c("black","black")) + theme_classic()

ggsave("dotplot_cluster_imp_gene_indiv_asp_norm_sortrebarplot1_half.svg", width=100, height=10, units="cm", dpi=96)
ggsave("dotplot_cluster_imp_gene_indiv_asp_norm_sortrebarplot1_half.jpg", width=100, height=10, units="cm", dpi=96)


#remove  Mrpl23-ps1
simprinted_gene_cords <- read.table("selected_imprinted_genes.cords.txt", header = F)
head(simprinted_gene_cords)
dim(simprinted_gene_cords)
simprinted_gene_cords <- simprinted_gene_cords[,c(1:3,5:6)]
colnames(simprinted_gene_cords) <- c("chr", "start", "end", "id", "Gene")
head(simprinted_gene_cords)
dim(simprinted_gene_cords)
head(gddsaSFNormcountsavg_chr_gene1_ratio)
dim(gddsaSFNormcountsavg_chr_gene1_ratio)
gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort = merge(gddsaSFNormcountsavg_chr_gene1_ratio, simprinted_gene_cords, by="id", all.x=FALSE)
head(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort)
dim(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort)
gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_bar <- gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort[,c(5,10,11)]
head(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_bar)
library(ggplot2)
library(gplots)
rownames(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_bar)
rownames(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_bar)=gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_bar[,1]
rownames(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_bar)
colnames(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_bar)
head(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_bar)
gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot = gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_bar[,-1]
gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot = as.matrix(t(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot))
head(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot)
dim(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot)
#barplot(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot, beside=TRUE, horiz=TRUE, col=c("red","blue","orange", "green"), xpd=FALSE,  xlim = c(0,1))


gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot1 <- stack(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot)
head(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot1)
colnames(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot1) <- c("Group", "Gene", "Ratio")
gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot1 <- data.frame(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot1)
head(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot1)
str(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot1)
gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot1_half <- gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot1
gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot1_half["Ratio"] <- gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot1_half$Ratio - 0.5
head(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot1_half)
library(ggplot2)
# Basic barplot
#way1 : ggplot2
pbar <-ggplot(data=gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot1_half, aes(x=Gene, y=Ratio, label=Group, color=Group, fill=Group)) +
  geom_bar(width = 0.5,stat="identity",  position=position_dodge(), size=0.6)
pbar + coord_flip()+
  scale_fill_manual(values=c("darkgreen","darkviolet"))+
  scale_color_manual(values=c("darkgreen","darkviolet")) + theme_classic()

#Apply fisher value 0.05, applying BH correction only on 36 imprinted genes leads to different fisheradjpvalue then applying on total genes. USE BH on total genes
gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher <- gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort
head(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher)
gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher["fisherpvalue"] <- apply(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher[,6:9],1, function(x) fisher.test(matrix(x,nrow =2))$p.value)
gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher["fisheradjpval"] <- p.adjust(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher$fisherpvalue,method="BH")
head(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher)
write.table(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher, "gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher.txt", sep="\t", quote = FALSE, append = FALSE)
gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05 <- gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher[which(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher$fisheradjpval <0.05),]
head(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05)
dim(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05)
write.table(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05, "gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05.txt", sep="\t", quote = FALSE, append = FALSE)

gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar <- gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05[,c(5,10,11)]
head(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar)
library(ggplot2)
library(gplots)
rownames(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar)
rownames(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar)=gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar[,1]
rownames(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar)
colnames(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar)
head(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar)
gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot = gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar[,-1]
gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot = as.matrix(t(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot))
head(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot)
dim(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot)
#barplot(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot, beside=TRUE, horiz=TRUE, col=c("red","blue","orange", "green"), xpd=FALSE,  xlim = c(0,1))
gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1 <- stack(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot)
head(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1)
colnames(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1) <- c("Group", "Gene", "Ratio")
gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1 <- data.frame(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1)
head(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1)
str(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1)
gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1_half <- gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1
gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1_half["Ratio"] <- gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1_half$Ratio - 0.5
head(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1_half)

library(ggplot2)
# Basic barplot
#way1 : ggplot2
pbar <-ggplot(data=gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1_half, aes(x=Gene, y=Ratio, label=Group, color=Group, fill=Group)) +
  geom_bar(width = 0.5,stat="identity",  position=position_dodge(), size=0.6)
pbar + coord_flip()+
  scale_fill_manual(values=c("darkgreen","darkviolet"))+
  scale_color_manual(values=c("darkgreen","darkviolet")) + theme_bw()
ggsave("gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1_half.svg", width=13, height=25, units="cm", dpi=96)
#theme(legend.position = "bottom", legend.direction = "horizontal", panel.background = element_blank(), axis.line.x = element_line(size = .6, colour = "black"))

#way2 : ggpubr
library(ggpubr)
ggbarplot(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_barplot1_half, x = "Gene", y = "Ratio", fill = "Group", color = "white", palette = c("darkgreen","darkviolet"),sort.by.groups = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Ratio",
          legend.title = "Expression",
          rotate = TRUE, position = position_dodge(),
          ggtheme = theme_classic())
ggbarplot(gddsaSFNormcountsavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1_half, x = "Gene", y = "Ratio", fill = "Group", color = "white", palette = c("darkgreen","darkviolet"),sort.by.groups = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Ratio",
          legend.title = "Expression",
          rotate = TRUE, position = position_dodge(),
          ggtheme = theme_classic())


######----------------------#######
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb")
data_imp_gene_indiv_asp_norm_part2 <- read.table("data_imp_gene_indiv_asp_norm_part2.txt", header = T)

# Basic barplot
head(data_imp_gene_indiv_asp_norm_part2)
dim(data_imp_gene_indiv_asp_norm_part2)
data_imp_gene_indiv_asp_norm_part2_are <- data_imp_gene_indiv_asp_norm_part2[,c(4,10:11,17:18)]
head(data_imp_gene_indiv_asp_norm_part2_are)
data_imp_gene_indiv_asp_norm_part2_sortre <- data_imp_gene_indiv_asp_norm_part2_are
head(data_imp_gene_indiv_asp_norm_part2_sortre)
rownames(data_imp_gene_indiv_asp_norm_part2_sortre) <- data_imp_gene_indiv_asp_norm_part2_sortre[,1]
data_imp_gene_indiv_asp_norm_part2_sortrebarplot <- data_imp_gene_indiv_asp_norm_part2_sortre[,-1]
data_imp_gene_indiv_asp_norm_part2_sortrebarplot <- as.matrix(t(data_imp_gene_indiv_asp_norm_part2_sortrebarplot))
head(data_imp_gene_indiv_asp_norm_part2_sortrebarplot)
dim(data_imp_gene_indiv_asp_norm_part2_sortrebarplot)

data_imp_gene_indiv_asp_norm_part2_sortrebarplot1 <- stack(data_imp_gene_indiv_asp_norm_part2_sortrebarplot)
head(data_imp_gene_indiv_asp_norm_part2_sortrebarplot1)
colnames(data_imp_gene_indiv_asp_norm_part2_sortrebarplot1) <- c("Group", "Gene", "Ratio")
data_imp_gene_indiv_asp_norm_part2_sortrebarplot1 <- data.frame(data_imp_gene_indiv_asp_norm_part2_sortrebarplot1)
head(data_imp_gene_indiv_asp_norm_part2_sortrebarplot1)
str(data_imp_gene_indiv_asp_norm_part2_sortrebarplot1)
data_imp_gene_indiv_asp_norm_part2_sortrebarplot1_half <- data_imp_gene_indiv_asp_norm_part2_sortrebarplot1
data_imp_gene_indiv_asp_norm_part2_sortrebarplot1_half["Ratio"] <- as.numeric(data_imp_gene_indiv_asp_norm_part2_sortrebarplot1_half$Ratio) - 0.5
head(data_imp_gene_indiv_asp_norm_part2_sortrebarplot1_half)
library(ggplot2)
# Basic dotplot

pdot2 <-ggplot(data=data_imp_gene_indiv_asp_norm_part2_sortrebarplot1_half, aes(x=Gene, y=Ratio, label=Group, color="black", fill=Group)) + 
  geom_hline(yintercept = c(-0.17,0,0.17), colour = "grey", linetype=c("dashed","solid","dashed")) + 
  geom_dotplot(binaxis='y', stackdir='centerwhole', stackratio=1.5, dotsize=1.2) + theme_classic() 
pdot2 + coord_flip()+ ylim(c(-0.5,0.5))+
  scale_fill_manual(values=c("green","green","red","red"))+
  scale_color_manual(values=c("black","black","black","black")) 

ggsave("dotplot_cluster_imp_gene_indiv_asp_norm_sortrebarplot2_half.svg", width=20, height=10, units="cm", dpi=96)
ggsave("dotplot_cluster_imp_gene_indiv_asp_norm_sortrebarplot2_half.jpg", width=20, height=10, units="cm", dpi=96)

cluster_imp_gene_indiv_asp_norm_part2 <- read.csv("cluster_imp_gene_indiv_asp_norm_part2.csv", header = T)
head(cluster_imp_gene_indiv_asp_norm_part2)
dim(cluster_imp_gene_indiv_asp_norm_part2)
summary(cluster_imp_gene_indiv_asp_norm_part2[,c(1,10:11,17:18)])

#Separate replicates
cluster_imp_gene_indiv_asp_norm_part2_are <- cluster_imp_gene_indiv_asp_norm_part2[,c(1,10:11,17:18)]
head(cluster_imp_gene_indiv_asp_norm_part2_are)
cluster_imp_gene_indiv_asp_norm_part2_sortre <- cluster_imp_gene_indiv_asp_norm_part2_are
head(cluster_imp_gene_indiv_asp_norm_part2_sortre)
rownames(cluster_imp_gene_indiv_asp_norm_part2_sortre) <- cluster_imp_gene_indiv_asp_norm_part2_sortre[,1]
cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot <- cluster_imp_gene_indiv_asp_norm_part2_sortre[,-1]
cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot <- as.matrix(t(cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot))
head(cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot)
dim(cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot)
cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1 <- stack(cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot)
head(cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1)
colnames(cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1) <- c("Group", "Gene", "Ratio")
cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1 <- data.frame(cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1)
head(cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1)
dim(cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1)
cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1["Group.1"] <- data.frame(rep(c(rep("WT_allele_ratio",2),rep("ZFP57KO_allele_ratio",2)), 2))
head(cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1)
dim(cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1)
str(cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1)
cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1_half <- cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1
cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1_half["Ratio"] <- as.numeric(cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1_half$Ratio) - 0.5
head(cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1_half)
library(ggplot2)
# Dotplot for Cluster
# Change the position
pcluster2 <-ggplot(cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1_half, aes(x=Gene, y=Ratio, fill= Group.1)) + geom_hline(yintercept = c(-0.17,0,0.17), colour = "grey", linetype=c("dashed","solid","dashed"))+
  geom_boxplot(position=position_dodge(0.8))+
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.8))
pcluster2 + scale_fill_manual(values=c("green","red"))+ylim(c(-0.5,0.5))+
  scale_color_manual(values=c("black","black")) + theme_classic()
#save manually
ggsave("dotplot_cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1_half.svg", width=20, height=10, units="cm", dpi=96)
ggsave("dotplot_cluster_imp_gene_indiv_asp_norm_part2_sortrebarplot1_half.jpg", width=20, height=10, units="cm", dpi=96)

#After Prof's suggestion only the statistically significant genes (fisher/BH adj p-value  < 0.1) were selected 
basilia_genes_for_dotplot <- read.table("basilia_genes_for_dotplot.txt")
colnames(basilia_genes_for_dotplot) <- c("Genes")
basilia_imp_genes_ddsSFNormcounts_chr.pos = merge(basilia_genes_for_dotplot, gddsaSFNormcounts_chr_gene1_fishercoordinate, by="Genes", all.x=FALSE)
head(basilia_imp_genes_ddsSFNormcounts_chr.pos)
dim(basilia_imp_genes_ddsSFNormcounts_chr.pos)

#Select Fisher/BH adj p-value < 0.05
basilia_imp_genes_ddsSFNormcounts_chr_sig <- basilia_imp_genes_ddsSFNormcounts_chr.pos[which(basilia_imp_genes_ddsSFNormcounts_chr.pos$fisheradjpval < 0.05),]
dim(basilia_imp_genes_ddsSFNormcounts_chr_sig)
head(basilia_imp_genes_ddsSFNormcounts_chr_sig)
basilia_imp_genes_ddsSFNormcounts_chr_resig <- basilia_imp_genes_ddsSFNormcounts_chr_sig
basilia_imp_genes_ddsSFNormcounts_chr_resig["WT_allele_ratio_1"] <- basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_WT_Rep1.JF1 / (basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_WT_Rep1.B6 + basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_WT_Rep1.JF1)
basilia_imp_genes_ddsSFNormcounts_chr_resig["WT_allele_ratio_2"] <- basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_WT_Rep2.JF1 / (basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_WT_Rep2.B6 + basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_WT_Rep2.JF1)
basilia_imp_genes_ddsSFNormcounts_chr_resig["WT_allele_ratio"] <- (basilia_imp_genes_ddsSFNormcounts_chr_resig$WT_allele_ratio_1 + basilia_imp_genes_ddsSFNormcounts_chr_resig$WT_allele_ratio_2)/2
basilia_imp_genes_ddsSFNormcounts_chr_resig["ZFP57KO_allele_ratio_1"] <- basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_ZFP57_KO_Rep1.JF1 / (basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_ZFP57_KO_Rep1.B6 + basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_ZFP57_KO_Rep1.JF1)
basilia_imp_genes_ddsSFNormcounts_chr_resig["ZFP57KO_allele_ratio_2"] <- basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_ZFP57_KO_Rep2.JF1 / (basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_ZFP57_KO_Rep2.B6 + basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_ZFP57_KO_Rep2.JF1)
basilia_imp_genes_ddsSFNormcounts_chr_resig["ZFP57KO_allele_ratio"] <- (basilia_imp_genes_ddsSFNormcounts_chr_resig$ZFP57KO_allele_ratio_1 + basilia_imp_genes_ddsSFNormcounts_chr_resig$ZFP57KO_allele_ratio_2)/2

data_imp_gene_indiv_asp_norm_resig <- basilia_imp_genes_ddsSFNormcounts_chr_resig
# order
data_imp_gene_indiv_asp_norm_resig <- data_imp_gene_indiv_asp_norm_resig[order(data_imp_gene_indiv_asp_norm_resig$WT_allele_ratio),]
head(data_imp_gene_indiv_asp_norm_resig)
dim(data_imp_gene_indiv_asp_norm_resig)
write.table(data_imp_gene_indiv_asp_norm_resig, "data_imp_gene_indiv_asp_norm_resig.txt", sep = "\t", quote = F, append = F, row.names = F)
head(data_imp_gene_indiv_asp_norm_resig,1)
dim(data_imp_gene_indiv_asp_norm_resig)
summary(data_imp_gene_indiv_asp_norm_resig[,c(20,21,23,24)])
data_imp_gene_indiv_asp_norm_resig_re <- data_imp_gene_indiv_asp_norm_resig[,c(1,20,21,23,24)]
head(data_imp_gene_indiv_asp_norm_resig_re)
rownames(data_imp_gene_indiv_asp_norm_resig_re) <- data_imp_gene_indiv_asp_norm_resig_re[,1]
data_imp_gene_indiv_asp_norm_resig_rebarplot <- data_imp_gene_indiv_asp_norm_resig_re[,-1]
data_imp_gene_indiv_asp_norm_resig_rebarplot <- as.matrix(t(data_imp_gene_indiv_asp_norm_resig_rebarplot))
head(data_imp_gene_indiv_asp_norm_resig_rebarplot)
dim(data_imp_gene_indiv_asp_norm_resig_rebarplot)
#barplot(data_imp_gene_indiv_asp_norm_resig_rebarplot, beside=TRUE, horiz=TRUE, col=c("red","blue","orange", "green"), xpd=FALSE,  xlim = c(0,1))
data_imp_gene_indiv_asp_norm_resig_rebarplot1 <- stack(data_imp_gene_indiv_asp_norm_resig_rebarplot)
head(data_imp_gene_indiv_asp_norm_resig_rebarplot1)

colnames(data_imp_gene_indiv_asp_norm_resig_rebarplot1) <- c("Group", "Gene", "Ratio")
data_imp_gene_indiv_asp_norm_resig_rebarplot1 <- data.frame(data_imp_gene_indiv_asp_norm_resig_rebarplot1)
head(data_imp_gene_indiv_asp_norm_resig_rebarplot1)
str(data_imp_gene_indiv_asp_norm_resig_rebarplot1)
data_imp_gene_indiv_asp_norm_resig_rebarplot1_half <- data_imp_gene_indiv_asp_norm_resig_rebarplot1
data_imp_gene_indiv_asp_norm_resig_rebarplot1_half["Ratio"] <- as.numeric(data_imp_gene_indiv_asp_norm_resig_rebarplot1_half$Ratio) - 0.5
head(data_imp_gene_indiv_asp_norm_resig_rebarplot1_half)
write.table(data_imp_gene_indiv_asp_norm_resig_rebarplot1_half, "data_imp_gene_indiv_asp_norm_resig_rebarplot1_half.txt", sep = "\t", quote = F, append = F)

psigdotsquare <-ggplot(data=data_imp_gene_indiv_asp_norm_resig_rebarplot1_half, aes(x=Gene, y=Ratio, label=Group, color="black", fill=Group)) + 
  geom_hline(yintercept = c(-0.17,0,0.17), colour = "grey", linetype=c("dashed","solid","dashed")) + 
  geom_point(aes(shape=Group, color=Group, fill=Group, size=Group), position = position_dodge(width = 0.1)) + theme_classic()  

psigdotsquare + coord_flip()+ ylim(c(-0.5,0.5))+
  scale_fill_manual(values=c("green","green","red","red"))+
  scale_color_manual(values=c("black","black","black","black")) +
  scale_shape_manual(values=c(22,24,22,24))+
  scale_size_manual(values=c(1.5,1.5,1.5,1.5)) 

ggsave("trisquareplot_data_imp_gene_indiv_asp_norm_resig_sortrebarplot1_half.svg", width=13, height=22, units="cm", dpi=96)
ggsave("trisquareplot_data_imp_gene_indiv_asp_norm_resig_sortrebarplot1_half.jpg", width=13, height=22, units="cm", dpi=96)



###################  End of Analysis ##################################à





############## CPM normalization ###############  cpm function #see last
#Counts per million normalization (non-filtered normalization)
colSums(Bulkcountdata) # 7663520  10225378  9545119  10337082
head(gcountdata)
gcountdata <- data.frame(gcountdata)
gcountdataCPM <- data.frame(cbind(data.frame(rownames(gcountdata)),
                                  (gcountdata$JB1_WT_Rep1.B6 *1000000)/7663520,
                                  (gcountdata$JB1_WT_Rep2.B6*1000000)/10225378,
                                  (gcountdata$JB1_WT_Rep1.JF1*1000000)/7663520,
                                  (gcountdata$JB1_WT_Rep2.JF1*1000000)/10225378,
                                  (gcountdata$JB1_ZFP57_KO_Rep1.B6*1000000)/9545119,
                                  (gcountdata$JB1_ZFP57_KO_Rep2.B6*1000000)/10337082,
                                  (gcountdata$JB1_ZFP57_KO_Rep1.JF1*1000000)/9545119,
                                  (gcountdata$JB1_ZFP57_KO_Rep2.JF1*1000000)/10337082), stringsAsFactors = FALSE)

head(gcountdataCPM)
dim(gcountdataCPM)
colnames(gcountdataCPM) <-  c("id","JB1_WT_Rep1.B6","JB1_WT_Rep2.B6","JB1_WT_Rep1.JF1","JB1_WT_Rep2.JF1","JB1_ZFP57_KO_Rep1.B6","JB1_ZFP57_KO_Rep2.B6","JB1_ZFP57_KO_Rep1.JF1","JB1_ZFP57_KO_Rep2.JF1")
rownames(gcountdataCPM) <- gcountdataCPM[,1]
gcountdataCPM <- gcountdataCPM[,-1]
head(gcountdataCPM)
gcountdataCPM <- data.frame(gcountdataCPM)
boxplot(gcountdataCPM[rowSums(gcountdataCPM) >= 10,], ylim=c(0,50))
plotPCA(as.matrix(gcountdataCPM), labels=F, col =  c("blue","blue","darkred","darkred","darkgreen","darkgreen","orange","orange"))
#PCA
gcountdataCPMfilt <- gcountdataCPM[rowSums(gcountdataCPM) >= 10,]
tgcountdataCPMfilt = t(gcountdataCPMfilt)
dim(tgcountdataCPMfilt)
rownames(tgcountdataCPMfilt)
tgcountdataCPMfilt = data.frame(tgcountdataCPMfilt)
tgcountdataCPMfilt["Color"] <-   c("WTB6_1","WTB6_2","WTJF1_1","WTJF1_2","ZFP57KOB6_1","ZFP57KOB6_2","ZFP57KOJF1_1","ZFP57KOJF1_2")
dim(tgcountdataCPMfilt)

dfx <-tgcountdataCPMfilt[c(1:10450)]
PC<-prcomp(dfx, scale. = T)
PCi<-data.frame(PC$x,Color=tgcountdataCPMfilt$Color)
percentage <- round(PC$sdev^2 / sum(PC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentage <- paste( colnames(PCi), "(", paste( as.character(percentage), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

p14<-ggplot(PCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentage[1]) + ylab(percentage[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("blue","blue","darkred","darkred","darkgreen","darkgreen","orange","orange"))+
  scale_shape_manual(values=c(19,17,19,17,19,17,19,17))
p14 <- p14+theme_bw()
p14

ggsave("PCA_tgcountdataCPMfilt_scaleT.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T

head(gcountdataCPM)
#Normalization in allele specific: https://support.bioconductor.org/p/73239/
gcountdataCPM_sorted = data.frame(gcountdataCPM[order(rownames(gcountdataCPM)),])
head(gcountdataCPM_sorted)
dim(gcountdataCPM_sorted)
gcountdataCPM_sorted.id <- cbind.data.frame(data.frame(rownames(gcountdataCPM_sorted)), gcountdataCPM_sorted)
head(gcountdataCPM_sorted.id) 
colnames(gcountdataCPM_sorted.id) <- c("id", "JB1_WT_Rep1.B6", "JB1_WT_Rep2.B6", "JB1_WT_Rep1.JF1", "JB1_WT_Rep2.JF1", "JB1_ZFP57_KO_Rep1.B6", "JB1_ZFP57_KO_Rep2.B6", "JB1_ZFP57_KO_Rep1.JF1", "JB1_ZFP57_KO_Rep2.JF1")
gcountdataCPM_sortedavg <- data.frame(cbind(data.frame(rownames(gcountdataCPM_sorted)),
                                            (gcountdataCPM_sorted$JB1_WT_Rep1.B6 + gcountdataCPM_sorted$JB1_WT_Rep2.B6)/2, 
                                            (gcountdataCPM_sorted$JB1_WT_Rep1.JF1 + gcountdataCPM_sorted$JB1_WT_Rep2.JF1)/2,
                                            (gcountdataCPM_sorted$JB1_ZFP57_KO_Rep1.B6 + gcountdataCPM_sorted$JB1_ZFP57_KO_Rep2.B6)/2,
                                            (gcountdataCPM_sorted$JB1_ZFP57_KO_Rep1.JF1 + gcountdataCPM_sorted$JB1_ZFP57_KO_Rep2.JF1)/2))

head(gcountdataCPM_sortedavg)
dim(gcountdataCPM_sortedavg)
colnames(gcountdataCPM_sortedavg) <-  c("id","JB1_WTB6","JB1_WTJF1","JB1_ZFP57KOB6","JB1_ZFP57KOJF1")
rownames(gcountdataCPM_sortedavg) <- gcountdataCPM_sortedavg[,1]
head(gcountdataCPM_sortedavg)
gcountdataCPM_sortedavg.id <- data.frame(gcountdataCPM_sortedavg$id)
head(gcountdataCPM_sortedavg.id)
colnames(gcountdataCPM_sortedavg.id) <- "id"
#z-score
head(gcountdataCPM_sortedavg)
dim(gcountdataCPM_sortedavg)
#Transform and scale
z_TgcountdataCPM_sortedavg= scale(t(gcountdataCPM_sortedavg[,2:5]), center = TRUE, scale = TRUE)
z_gcountdataCPM_sortedavg <- t(z_TgcountdataCPM_sortedavg)
head(z_gcountdataCPM_sortedavg)
dim(z_gcountdataCPM_sortedavg)
write.table(z_gcountdataCPM_sortedavg, "z_gcountdataCPM_sortedavg.txt", sep="\t", quote=F, col.names=NA)

z_gcountdataCPM_sortedavg.id <- data.frame(cbind(gcountdataCPM_sortedavg.id, z_gcountdataCPM_sortedavg))
head(z_gcountdataCPM_sortedavg.id)
tail(z_gcountdataCPM_sortedavg.id)
write.table(z_gcountdataCPM_sortedavg.id, "z_gcountdataCPM_sortedavg.id.txt", sep="\t", quote=F, col.names=T, row.names = F)

chr.pos.CPMavg =  read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt",header=FALSE)
head(chr.pos.CPMavg)
colnames(chr.pos.CPMavg) <- c("id", "Genes", "chr", "start", "end")
head(chr.pos.CPMavg)
chr.pos.CPMavg.2 = merge(gcountdataCPM_sortedavg.id, chr.pos.CPMavg, by="id", all.x=TRUE)
head(chr.pos.CPMavg.2)
dim(chr.pos.CPMavg.2)
gcountdataCPM_sortedavg_chr_gene <- data.frame(cbind(chr.pos.CPMavg.2, gcountdataCPM_sortedavg))
head(gcountdataCPM_sortedavg_chr_gene)
gcountdataCPM_sortedavg_chr_gene1 <- gcountdataCPM_sortedavg_chr_gene[,c(3:5,2,1,6:10)]
head(gcountdataCPM_sortedavg_chr_gene1)
write.table(gcountdataCPM_sortedavg_chr_gene1, "gcountdataCPM_sortedavg_chr_gene1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)

chr.pos.CPM =  read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt",header=FALSE)
head(chr.pos.CPM)
colnames(chr.pos.CPM) <- c("id", "Genes", "chr", "start", "end")
head(chr.pos.CPM)
head(gcountdataCPM_sorted.id)
chr.pos.CPM.2 = merge(gcountdataCPM_sorted.id, chr.pos.CPM, by="id", all.x=TRUE)
head(chr.pos.CPM.2)
dim(chr.pos.CPM.2)
head(gcountdataCPM_sorted)
gcountdataCPM_sorted_chr_gene <- data.frame(cbind(chr.pos.CPM.2, gcountdataCPM_sorted))
head(gcountdataCPM_sorted_chr_gene)
gcountdataCPM_sorted_chr_gene1 <- gcountdataCPM_sorted_chr_gene[,c(11:13,10,1,2:9)]
head(gcountdataCPM_sorted_chr_gene1)
write.table(gcountdataCPM_sorted_chr_gene1, "gcountdataCPM_sorted_chr_gene1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)
dim(gcountdataCPM_sorted_chr_gene1)

chr.pos.CPM.3 = merge(z_gcountdataCPM_sortedavg.id, chr.pos.CPM, by="id", all.x=TRUE)
head(chr.pos.CPM.3)
dim(chr.pos.CPM.3)
z_gcountdataCPM_sortedavg_chr_gene <- data.frame(cbind(chr.pos.CPM.3, z_gcountdataCPM_sortedavg))
head(z_gcountdataCPM_sortedavg_chr_gene)
z_gcountdataCPM_sortedavg_chr_gene1 <- z_gcountdataCPM_sortedavg_chr_gene[,c(7:9,6,1,2:5)]
head(z_gcountdataCPM_sortedavg_chr_gene1)
tail(z_gcountdataCPM_sortedavg_chr_gene1)
dim(z_gcountdataCPM_sortedavg_chr_gene1)
write.table(z_gcountdataCPM_sortedavg_chr_gene1, "z_gcountdataCPM_sortedavg_chr_gene1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)


#Individuals replicates separate process
head(gcountdataCPM)
gcountdataCPM["id"] <- data.frame(rownames(gcountdataCPM))
head(gcountdataCPM)
chr.pos.CPM.4 = merge(gcountdataCPM, chr.pos.CPM, by="id", all.x=TRUE)
head(chr.pos.CPM.4)
dim(chr.pos.CPM.4)
gcountdataCPM_chr_gene <- chr.pos.CPM.4
head(gcountdataCPM_chr_gene)
gcountdataCPM_chr_gene1 <- gcountdataCPM_chr_gene[,c(11:13,10,1,2:9)]
head(gcountdataCPM_chr_gene1)
write.table(gcountdataCPM_chr_gene1, "gcountdataCPM_chr_gene1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)

#Get imprinted genes
#This file imprinted_gene_name.dups.txt was produced from my search and bouschet imprinted genes and duplicates were removed
sort -k1,1 imprinted_gene_name.dups.txt -u > imprinted_gene_name.txt #188
/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt

#Extract imprinted genes data
#Selected genes list was given by Basilia (36) and coordinates were obtained by fgrep -f with gtf gene file. 
#Some genes like Inpp5f and Ddc were filtere out because they were not found to be imprinted leaving 34 genes 
#I removed it before itself as they can affect the z-score scaling so I kept only 34 genes
selectedgenes.cords <- read.table("selected_imprinted_genes_refilt.cords.txt", header = F)
head(selectedgenes.cords)
selectedgenes.cords<- selectedgenes.cords[,c(1:3,5:6)]
colnames(selectedgenes.cords) <- c("chr", "start", "end", "id", "Gene")
z_gcountdataCPM_sortedavg.selected_genes = merge(z_gcountdataCPM_sortedavg.id, selectedgenes.cords, by="id", all.x=FALSE)
head(z_gcountdataCPM_sortedavg.selected_genes)
dim(z_gcountdataCPM_sortedavg.selected_genes)
library(ggplot2)
library(gplots)
#setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb")
#data <- read.table("all_DEG_rearranged_filt_ImpGenes.heatmap.txt",header=TRUE)
z_gcountdataCPM_sortedavg.selected_genes_sort <- z_gcountdataCPM_sortedavg.selected_genes[order(z_gcountdataCPM_sortedavg.selected_genes$chr, z_gcountdataCPM_sortedavg.selected_genes$start),]
head(z_gcountdataCPM_sortedavg.selected_genes_sort)
dim(z_gcountdataCPM_sortedavg.selected_genes_sort)
z_gcountdataCPM_sortedavg.selected_genes_sort <- z_gcountdataCPM_sortedavg.selected_genes_sort[c(1,13:34, 2:12),]
head(z_gcountdataCPM_sortedavg.selected_genes_sort)
z_gcountdataCPM_sortedavg.selected_genes_sort_heat <- z_gcountdataCPM_sortedavg.selected_genes_sort[,c(9,2:5)]
rownames(z_gcountdataCPM_sortedavg.selected_genes_sort_heat)
rownames(z_gcountdataCPM_sortedavg.selected_genes_sort_heat)=z_gcountdataCPM_sortedavg.selected_genes_sort_heat[,1]
rownames(z_gcountdataCPM_sortedavg.selected_genes_sort_heat)
colnames(z_gcountdataCPM_sortedavg.selected_genes_sort_heat)
head(z_gcountdataCPM_sortedavg.selected_genes_sort_heat)
z_gcountdataCPM_sortedavg.selected_genes_sort_heat1 = z_gcountdataCPM_sortedavg.selected_genes_sort_heat[,-1]
z_gcountdataCPM_sortedavg.selected_genes_sort_heat1 = as.matrix(z_gcountdataCPM_sortedavg.selected_genes_sort_heat1)
head(z_gcountdataCPM_sortedavg.selected_genes_sort_heat1)
dim(z_gcountdataCPM_sortedavg.selected_genes_sort_heat1)
write.table(z_gcountdataCPM_sortedavg.selected_genes_sort_heat1, "z_gcountdataCPM_sortedavg.selected_genes_sort_heat1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)

#colfunc <- colorRampPalette( colors = brewer.pal(9,"PRGn") )
colfunc <- colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","white","#FDDBC7","#F4A582","#D6604D","#B2182B"))

#heatmap.2(data1, Colv = "NA",trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=1, font=3, cexCol = 1, margins =c(3,8), srtCol = 45, breaks = seq(-5,5, length.out = 100))
svg(filename="z_gcountdataCPM_sortedavg.selected_genes_sort_heat1.svg", width=8, height=18, pointsize=12)
heatmap.2(z_gcountdataCPM_sortedavg.selected_genes_sort_heat1, Colv = "NA",Rowv ='NA',trace = "none", col = colfunc ,
          lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5),
          density.info=c("none"), dendrogram="none", scale = "none", 
          sepwidth=c(0.001, 0.001), cexRow=2, font=3, cexCol = 0.8, 
          sepcolor="black", margins =c(6,8), srtCol = 45, 
          breaks = seq(-1,1, length.out = 100), colsep=1:ncol(z_gcountdataCPM_sortedavg.selected_genes_sort_heat1),
          rowsep=1:nrow(z_gcountdataCPM_sortedavg.selected_genes_sort_heat1))
dev.off()

library(pheatmap)
library(RColorBrewer)
breaksList = seq(-1, 1)
#pheatmap(data1,treeheight_row = 0, cluster_cols=F, cluster_rows=F, treeheight_col = 0, gaps_col =NULL, gaps_row = NULL, border_color = "black", breaks = breaksList,color= colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)))
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_gcountdataCPM_sortedavg.selected_genes_sort_heat1,
         color = colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","white","#FDDBC7","#F4A582","#D6604D","#B2182B"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=F, 
         cluster_rows=F,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12, 
         filename = "z_gcountdataCPM_sortedavg.selected_genes_sort_heat1.png")

dev.off()
#Now filter the matrix to remove low expressed genes, so cpm counts >=rowsum 10 are removed
#write.table(results, "results.txt", sep="\t", quote=F, col.names=NA)
gcountdataCPM_sortedavg.id.avg <- data.frame(gcountdataCPM_sortedavg)
head(gcountdataCPM_sortedavg.id.avg)
tail(gcountdataCPM_sortedavg.id.avg)
gcountdataCPM_sortedavg.selected_genes = merge(gcountdataCPM_sortedavg.id.avg, selectedgenes.cords, by="id", all.x=FALSE)
head(gcountdataCPM_sortedavg.selected_genes)
gcountdataCPM_sortedavg.selected_genes_sort <- gcountdataCPM_sortedavg.selected_genes[order(gcountdataCPM_sortedavg.selected_genes$chr, gcountdataCPM_sortedavg.selected_genes$start),]
head(gcountdataCPM_sortedavg.selected_genes_sort)
dim(gcountdataCPM_sortedavg.selected_genes_sort)
gcountdataCPM_sortedavg.selected_genes_sort <- gcountdataCPM_sortedavg.selected_genes_sort[c(1,13:34, 2:12),]
head(gcountdataCPM_sortedavg.selected_genes_sort)
gcountdataCPM_sortedavg.selected_genes_sort_heat1 <- gcountdataCPM_sortedavg.selected_genes_sort[,c(9,2:5)]
rownames(gcountdataCPM_sortedavg.selected_genes_sort_heat1)
rownames(gcountdataCPM_sortedavg.selected_genes_sort_heat1)=gcountdataCPM_sortedavg.selected_genes_sort_heat1[,1]
rownames(gcountdataCPM_sortedavg.selected_genes_sort_heat1)
head(gcountdataCPM_sortedavg.selected_genes_sort_heat1)
write.table(gcountdataCPM_sortedavg.selected_genes_sort_heat1, "gcountdataCPM_sortedavg.selected_genes_sort_heat1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb")

gcountdataCPM_sortedavg.selected_genes_sort_heat1 <- read.table("gcountdataCPM_sortedavg.selected_genes_sort_heat1.txt", header = T)
head(gcountdataCPM_sortedavg.selected_genes_sort_heat1)
gcountdataCPM_sortedavg.selected_genes_sort_heat1["WTratio"] <- gcountdataCPM_sortedavg.selected_genes_sort_heat1$JB1_WTJF1 / (gcountdataCPM_sortedavg.selected_genes_sort_heat1$JB1_WTB6 + gcountdataCPM_sortedavg.selected_genes_sort_heat1$JB1_WTJF1)
gcountdataCPM_sortedavg.selected_genes_sort_heat1["ZFP57KOratio"]  <- gcountdataCPM_sortedavg.selected_genes_sort_heat1$JB1_ZFP57KOJF1 / (gcountdataCPM_sortedavg.selected_genes_sort_heat1$JB1_ZFP57KOB6 + gcountdataCPM_sortedavg.selected_genes_sort_heat1$JB1_ZFP57KOJF1)
gcountdataCPM_sortedavg.selected_genes_sort_heat1["B6ratio"] <- log2(gcountdataCPM_sortedavg.selected_genes_sort_heat1$JB1_ZFP57KOB6 / gcountdataCPM_sortedavg.selected_genes_sort_heat1$JB1_WTB6)
gcountdataCPM_sortedavg.selected_genes_sort_heat1["JF1ratio"]  <- log2(gcountdataCPM_sortedavg.selected_genes_sort_heat1$JB1_ZFP57KOJF1 / gcountdataCPM_sortedavg.selected_genes_sort_heat1$JB1_WTJF1)
write.table(gcountdataCPM_sortedavg.selected_genes_sort_heat1, "gcountdataCPM_sortedavg.selected_genes_sort_heat1.withratio.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)

head(gcountdataCPM_sortedavg.selected_genes_sort_heat1)
library(pheatmap)
library(RColorBrewer)
#Predict Mnoallecity or Biallecity
#pheatmap(data1,treeheight_row = 0, cluster_cols=F, cluster_rows=F, treeheight_col = 0, gaps_col =NULL, gaps_row = NULL, border_color = "black", breaks = breaksList,color= colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)))
breaksList1 = c(0, 0.31, 0.67)
pheatmap(gcountdataCPM_sortedavg.selected_genes_sort_heat1[,6:7],
         color = colorRampPalette(c("#00B0F0","#f8de38","#f838b0"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "black",
         cluster_cols=F, 
         cluster_rows=F,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12,
         filename = "z_gcountdataCPM_sortedavg.selected_genes_sort_heat1_Interpretation.png")
dev.off()
#Predict silent or active allele

#Copy gcountdataCPM_sortedavg.selected_genes_sort_heat1.withratio.txt file in excel and separate silent and active allele and save as gcountdataCPM_sortedavg.selected_genes_sort_heat1.withratio.reorder.txt. The excel sheet is in 5mb folder names as order heatmap silent active.xlsx
gcountdataCPM_sortedavg.selected_genes_sort_heat1.withratio.reorder <- read.table("gcountdataCPM_sortedavg.selected_genes_sort_heat1.withratio.reorder.txt", header = T)
head(gcountdataCPM_sortedavg.selected_genes_sort_heat1.withratio.reorder)
rownames(gcountdataCPM_sortedavg.selected_genes_sort_heat1.withratio.reorder) <- gcountdataCPM_sortedavg.selected_genes_sort_heat1.withratio.reorder[,1]
gcountdataCPM_sortedavg.selected_genes_sort_heat1.withratio.reorder <- gcountdataCPM_sortedavg.selected_genes_sort_heat1.withratio.reorder[,-1]
gcountdataCPM_sortedavg.selected_genes_sort_heat1.withratio.reorder
breaksList1 = c(-300,-200,-0.15,0.15)
pheatmap(gcountdataCPM_sortedavg.selected_genes_sort_heat1.withratio.reorder[,1:2],
         color = colorRampPalette(c("darkgrey","#a10517","white","#35b04b"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "black",
         cluster_cols=F, 
         cluster_rows=F,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12,
         filename = "z_gcountdataCPM_sortedavg.selected_genes_sort_heat1_Interpretation_updown.png")

dev.off()
#Heatmap for imprinted genes
grep ENS imprinted_gene_rev.sort.txt > imprinted_gene_cords.bed #grep coordinates

imprinted_gene_cords <- read.table("imprinted_gene_cords.bed", header = F)
head(imprinted_gene_cords)
colnames(imprinted_gene_cords) <- c("chr", "start", "end", "id", "Gene")
head(imprinted_gene_cords)
dim(imprinted_gene_cords)
head(z_gcountdataCPM_sortedavg)
dim(z_gcountdataCPM_sortedavg)
z_gcountdataCPM_sortedavg_sort <- data.frame(cbind(rownames(z_gcountdataCPM_sortedavg), z_gcountdataCPM_sortedavg))
head(z_gcountdataCPM_sortedavg_sort)
colnames(z_gcountdataCPM_sortedavg_sort) <- c("id", "JB1_WTB6", "JB1_WTJF1", "JB1_ZFP57KOB6", "JB1_ZFP57KOJF1")
z_gcountdataCPM_sortedavg_imprinted_sort = merge(z_gcountdataCPM_sortedavg_sort, imprinted_gene_cords, by="id", all.x=F)
head(z_gcountdataCPM_sortedavg_imprinted_sort)
dim(z_gcountdataCPM_sortedavg_imprinted_sort)
write.table(z_gcountdataCPM_sortedavg_imprinted_sort, "z_gcountdataCPM_sortedavg_imprinted_sort.txt", quote = F, append = F )
library(ggplot2)
library(gplots)
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb")
z_gcountdataCPM_sortedavg_imprinted_sortre <- z_gcountdataCPM_sortedavg_imprinted_sort[,c(6:9,1:5)]
head(z_gcountdataCPM_sortedavg_imprinted_sortre)
write.table(z_gcountdataCPM_sortedavg_imprinted_sortre, "z_gcountdataCPM_sortedavg_imprinted_sortre.txt", quote = F, append= F, row.names = F)
z_gCPM_avgimppos <- read.table("z_gcountdataCPM_sortedavg_imprinted_sortre.txt",header=T)
head(z_gCPM_avgimppos)
rownames(z_gCPM_avgimppos) <- z_gCPM_avgimppos$Gene
z_gCPM_avgimppos1 = z_gCPM_avgimppos[,6:9]
z_gCPM_avgimppos2 = as.matrix(z_gCPM_avgimppos1)
head(z_gCPM_avgimppos2)
dim(z_gCPM_avgimppos2)

colfunc <- colorRampPalette(c("navy","white", "#C71111"))
#heatmap.2(z_gCPM_avgimp1, Colv = "NA",trace = "none", col = colfunc , lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5), density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=1, font=3, cexCol = 1, margins =c(3,8), srtCol = 45, breaks = seq(-5,5, length.out = 100))
heatmap.2(z_gCPM_avgimppos2, Colv = "NA",Rowv ='NA',trace = "none", col = colfunc ,
          lmat=rbind( c(5, 4, 2), c(6,1,3 ) ), lwid=c(1, 3, 2 ),lhei = c(1,7), keysize=1.2, key.par = list(cex=0.5),
          density.info=c("none"), dendrogram="none", scale = "none", 
          sepwidth=c(0.001, 0.001), cexRow=0.4, font=3, cexCol = 0.5, 
          sepcolor="black", margins =c(6,8), srtCol = 45, 
          breaks = seq(-2,2, length.out = 100), colsep=1:ncol(z_gCPM_avgimppos2),
          rowsep=1:nrow(z_gCPM_avgimppos2))

#save as z_gCPM_avgimppos2.svg
# load package
library(pheatmap)
library(RColorBrewer)
breaksList = seq(-3, 3)


pheatmap(z_gCPM_avgimppos2,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)),
         breaks = breaksList,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 1,
         cellwidth = 12, 
         cellheight = 30,
         filename = "z_gCPM_avgimppos2.png")

#save as z_gCPM_avgimppos2.svg

#Fisher test and tag chromsome positions
head(gcountdataCPM_sortedavg_chr_gene1)
dim(gcountdataCPM_sortedavg_chr_gene1)
gcountdataCPM_sortedavg_chr_gene1_fisher <- gcountdataCPM_sortedavg_chr_gene1[,7:10]
head(gcountdataCPM_sortedavg_chr_gene1_fisher)
gcountdataCPM_sortedavg_chr_gene1_fisher1 <- apply(gcountdataCPM_sortedavg_chr_gene1_fisher,1, function(x) fisher.test(matrix(x,nrow =2))$p.value)
#warning will be observed because we need only integers for runnung this function.
#In fisher.test(matrix(x, nrow = 2)) :'x' has been rounded to integer: Mean relative difference:
head(gcountdataCPM_sortedavg_chr_gene1_fisher1)
gcountdataCPM_sortedavg_chr_gene1_fisher["fisherpvalue"] <- data.frame(gcountdataCPM_sortedavg_chr_gene1_fisher1)
head(gcountdataCPM_sortedavg_chr_gene1_fisher)
gcountdataCPM_sortedavg_chr_gene1_fisher["fisheradjpval"] <- p.adjust(gcountdataCPM_sortedavg_chr_gene1_fisher1,method="BH")
head(gcountdataCPM_sortedavg_chr_gene1_fisher)
write.table(gcountdataCPM_sortedavg_chr_gene1_fisher, "gcountdataCPM_sortedavg_chr_gene1_fisher.txt", sep="\t", quote = FALSE, append = FALSE)
gcountdataCPM_sortedavg_chr_gene1_fisher_0.05 <- gcountdataCPM_sortedavg_chr_gene1_fisher[which(gcountdataCPM_sortedavg_chr_gene1_fisher$fisheradjpval <0.05),]
head(gcountdataCPM_sortedavg_chr_gene1_fisher_0.05)
dim(gcountdataCPM_sortedavg_chr_gene1_fisher_0.05)
write.table(gcountdataCPM_sortedavg_chr_gene1_fisher_0.05, "gcountdataCPM_sortedavg_chr_gene1_fisher_0.05.txt", sep="\t", quote = FALSE, append = FALSE)

#All imprinted genes maternal ratio
head(gcountdataCPM_sortedavg_chr_gene1)
gcountdataCPM_sortedavg_chr_gene1_imprinted_sort = merge(gcountdataCPM_sortedavg_chr_gene1, imprinted_gene_cords, by="id", all.x=F)
head(gcountdataCPM_sortedavg_chr_gene1_imprinted_sort)
dim(gcountdataCPM_sortedavg_chr_gene1_imprinted_sort)
gcountdataCPM_sortedavg_chr_gene1_imprinted_sortre <- gcountdataCPM_sortedavg_chr_gene1_imprinted_sort[,c(2:4,1,5,7:10)]
head(gcountdataCPM_sortedavg_chr_gene1_imprinted_sortre)
gcountdataCPM_sortedavg_chr_gene1_imprinted_sortre["WTmatratio"] <- gcountdataCPM_sortedavg_chr_gene1_imprinted_sortre$JB1_WTJF1 / (gcountdataCPM_sortedavg_chr_gene1_imprinted_sortre$JB1_WTB6 + gcountdataCPM_sortedavg_chr_gene1_imprinted_sortre$JB1_WTJF1)
gcountdataCPM_sortedavg_chr_gene1_imprinted_sortre["Zfp57komatratio"] <- gcountdataCPM_sortedavg_chr_gene1_imprinted_sortre$JB1_ZFP57KOJF1 / (gcountdataCPM_sortedavg_chr_gene1_imprinted_sortre$JB1_ZFP57KOB6 + gcountdataCPM_sortedavg_chr_gene1_imprinted_sortre$JB1_ZFP57KOJF1)
head(gcountdataCPM_sortedavg_chr_gene1_imprinted_sortre)
write.table(gcountdataCPM_sortedavg_chr_gene1_imprinted_sortre, "gcountdataCPM_sortedavg_chr_gene1_imprinted_sortre.txt", quote = F, append = F, row.names = F)

#grep chrX gcountdataCPM_sortedavg_chr_gene1_imprinted_sortre.txt -v | awk '{if($10 <= 0.31 || $10 >= 0.67) print $0}' > gcountdataCPM_sortedavg_chr_gene1_imprinted_sortre.filt.txt

#Scatter Plot
head(gcountdataCPM_sortedavg_chr_gene1)
gcountdataCPM_sortedavg_chr_gene1_ratio <- gcountdataCPM_sortedavg_chr_gene1
gcountdataCPM_sortedavg_chr_gene1_ratio["matWTratio"] <- gcountdataCPM_sortedavg_chr_gene1_ratio$JB1_WTJF1/(gcountdataCPM_sortedavg_chr_gene1_ratio$JB1_WTB6 + gcountdataCPM_sortedavg_chr_gene1_ratio$JB1_WTJF1)
head(gcountdataCPM_sortedavg_chr_gene1_ratio)
gcountdataCPM_sortedavg_chr_gene1_ratio["matZFP57KOratio"] <- gcountdataCPM_sortedavg_chr_gene1_ratio$JB1_ZFP57KOJF1/(gcountdataCPM_sortedavg_chr_gene1_ratio$JB1_ZFP57KOB6 + gcountdataCPM_sortedavg_chr_gene1_ratio$JB1_ZFP57KOJF1)
head(gcountdataCPM_sortedavg_chr_gene1_ratio)

plot(gcountdataCPM_sortedavg_chr_gene1_ratio$matWTratio, col="darkgreen", pch =16, cex = 0.3)
par(new=T)
plot(gcountdataCPM_sortedavg_chr_gene1_ratio$matZFP57KOratio, col="darkred", pch =16, cex = 0.3)

gcountdataCPM_sortedavg_chr_gene1_ratioR <- data.frame(gcountdataCPM_sortedavg_chr_gene1_ratio[,10:11])
head(gcountdataCPM_sortedavg_chr_gene1_ratioR[c(4500:4900),])
gcountdataCPM_sortedavgratioRange <- gcountdataCPM_sortedavg_chr_gene1_ratioR[c(4500:4900),]
head(gcountdataCPM_sortedavgratioRange)
plot(gcountdataCPM_sortedavgratioRange$WTratio, col="darkgreen", pch =16, cex = 0.8)
par(new=T)
plot(gcountdataCPM_sortedavgratioRange$ZFP57KOratio, col="darkred", pch =16, cex = 0.8)



#Barplot of 35 imprinted genes.

#Basilia provided the an excel sheet file which I saved as Table_for_ankit _27_july.csv.
#For CPM also I prepared the similar sheet using gcountdataCPM_chr_gene1.txt.
#fgrep -f basilia_genes_for_dotplot.txt gcountdataCPM_chr_gene1.txt -w | sort -k4,4 > basilia_genes_for_dotplot_gcountdataCPM.txt

#Manually done
#Open basilia_genes_for_dotplot_gcountdataCPM.txt and if 0 occurs in both allele, we replaced it with 1 (atleast one read), the other allele was still high
#save as basilia_genes_for_dotplot_gcountdataCPM_re.txt
basilia_genes_for_dotplot_gcountdataCPM <- read.table("basilia_genes_for_dotplot_gcountdataCPM_re.txt", header = F)

colnames(basilia_genes_for_dotplot_gcountdataCPM) <- colnames(gcountdataCPM_chr_gene1)
basilia_genes_for_dotplot_gcountdataCPM["WT_allele_ratio_1"] <- basilia_genes_for_dotplot_gcountdataCPM$JB1_WT_Rep1.JF1 / (basilia_genes_for_dotplot_gcountdataCPM$JB1_WT_Rep1.B6 + basilia_genes_for_dotplot_gcountdataCPM$JB1_WT_Rep1.JF1)
basilia_genes_for_dotplot_gcountdataCPM["WT_allele_ratio_2"] <- basilia_genes_for_dotplot_gcountdataCPM$JB1_WT_Rep2.JF1 / (basilia_genes_for_dotplot_gcountdataCPM$JB1_WT_Rep2.B6 + basilia_genes_for_dotplot_gcountdataCPM$JB1_WT_Rep2.JF1)
basilia_genes_for_dotplot_gcountdataCPM["WT_allele_ratio"] <- (basilia_genes_for_dotplot_gcountdataCPM$WT_allele_ratio_1 + basilia_genes_for_dotplot_gcountdataCPM$WT_allele_ratio_2)/2
basilia_genes_for_dotplot_gcountdataCPM["ZFP57KO_allele_ratio_1"] <- basilia_genes_for_dotplot_gcountdataCPM$JB1_ZFP57_KO_Rep1.JF1 / (basilia_genes_for_dotplot_gcountdataCPM$JB1_ZFP57_KO_Rep1.B6 + basilia_genes_for_dotplot_gcountdataCPM$JB1_ZFP57_KO_Rep1.JF1)
basilia_genes_for_dotplot_gcountdataCPM["ZFP57KO_allele_ratio_2"] <- basilia_genes_for_dotplot_gcountdataCPM$JB1_ZFP57_KO_Rep2.JF1 / (basilia_genes_for_dotplot_gcountdataCPM$JB1_ZFP57_KO_Rep2.B6 + basilia_genes_for_dotplot_gcountdataCPM$JB1_ZFP57_KO_Rep2.JF1)
basilia_genes_for_dotplot_gcountdataCPM["ZFP57KO_allele_ratio"] <- (basilia_genes_for_dotplot_gcountdataCPM$ZFP57KO_allele_ratio_1 + basilia_genes_for_dotplot_gcountdataCPM$ZFP57KO_allele_ratio_2)/2

gCPM_imp_gene_indiv_asp_norm <- basilia_genes_for_dotplot_gcountdataCPM
# order
gCPM_imp_gene_indiv_asp_norm <- gCPM_imp_gene_indiv_asp_norm[order(gCPM_imp_gene_indiv_asp_norm$WT_allele_ratio),]

head(gCPM_imp_gene_indiv_asp_norm)
dim(gCPM_imp_gene_indiv_asp_norm)
summary(gCPM_imp_gene_indiv_asp_norm[,c(16,19)])
gCPM_imp_gene_indiv_asp_norm_re <- gCPM_imp_gene_indiv_asp_norm[,c(4,16,19)]
head(gCPM_imp_gene_indiv_asp_norm_re)

colnames(gCPM_imp_gene_indiv_asp_norm_re) <- c("Gene", "WT_avg_allele_ratio", "ZFP57KO_avg_allele_ratio")
rownames(gCPM_imp_gene_indiv_asp_norm_re) <- gCPM_imp_gene_indiv_asp_norm_re[,1]
gCPM_imp_gene_indiv_asp_norm_rebarplot <- gCPM_imp_gene_indiv_asp_norm_re[,-1]
gCPM_imp_gene_indiv_asp_norm_rebarplot <- as.matrix(t(gCPM_imp_gene_indiv_asp_norm_rebarplot))
head(gCPM_imp_gene_indiv_asp_norm_rebarplot)
dim(gCPM_imp_gene_indiv_asp_norm_rebarplot)
#barplot(gCPM_imp_gene_indiv_asp_norm_rebarplot, beside=TRUE, horiz=TRUE, col=c("red","blue","orange", "green"), xpd=FALSE,  xlim = c(0,1))

gCPM_imp_gene_indiv_asp_norm_rebarplot1 <- stack(gCPM_imp_gene_indiv_asp_norm_rebarplot)
head(gCPM_imp_gene_indiv_asp_norm_rebarplot1)
colnames(gCPM_imp_gene_indiv_asp_norm_rebarplot1) <- c("Group", "Gene", "Ratio")
gCPM_imp_gene_indiv_asp_norm_rebarplot1 <- data.frame(gCPM_imp_gene_indiv_asp_norm_rebarplot1)
head(gCPM_imp_gene_indiv_asp_norm_rebarplot1)
str(gCPM_imp_gene_indiv_asp_norm_rebarplot1)
gCPM_imp_gene_indiv_asp_norm_rebarplot1_half <- gCPM_imp_gene_indiv_asp_norm_rebarplot1
gCPM_imp_gene_indiv_asp_norm_rebarplot1_half["Ratio"] <- as.numeric(gCPM_imp_gene_indiv_asp_norm_rebarplot1_half$Ratio) - 0.5
head(gCPM_imp_gene_indiv_asp_norm_rebarplot1_half)
write.table(gCPM_imp_gene_indiv_asp_norm_rebarplot1_half, "gCPM_imp_gene_indiv_asp_norm_rebarplot1_half.txt", sep = "\t", quote = F, append = F)
library(ggplot2)
# Basic barplot

#way1 : ggplot2
pbar <-ggplot(data=gCPM_imp_gene_indiv_asp_norm_rebarplot1_half, aes(x=Gene, y=Ratio, label=Group, color=Group, fill=Group)) +
  geom_bar(width = 0.5,stat="identity",  position=position_dodge(), size=0.6)
pbar + coord_flip()+
  scale_fill_manual(values=c("darkgreen","darkred"))+
  scale_color_manual(values=c("darkgreen","darkred")) + theme_bw()
ggsave("barplot_gCPM_imp_gene_indiv_asp_norm_rebarplot1_half.svg", width=13, height=25, units="cm", dpi=96)
ggsave("barplot_gCPM_imp_gene_indiv_asp_norm_rebarplot1_half.jpg", width=13, height=25, units="cm", dpi=96)


#way2 : ggpubr
library(ggpubr)
ggbarplot(gCPM_imp_gene_indiv_asp_norm_rebarplot1_half, x = "Gene", y = "Ratio", fill = "Group", color = "white", palette = c("darkgreen","darkred"),sort.by.groups = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Ratio",
          legend.title = "Expression",
          rotate = TRUE, position = position_dodge(),
          ggtheme = theme_classic())
ggbarplot(gCPM_imp_gene_indiv_asp_norm_rebarplot1_half, x = "Gene", y = "Ratio", fill = "Group", color = "white", palette = c("darkgreen","darkred"),sort.by.groups = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Ratio",
          legend.title = "Expression",
          rotate = TRUE, position = position_dodge(),
          ggtheme = theme_classic())

head(gCPM_imp_gene_indiv_asp_norm)
gCPM_imp_gene_indiv_asp_norm_are <- gCPM_imp_gene_indiv_asp_norm[,c(4,14:15,17:18)]
head(gCPM_imp_gene_indiv_asp_norm_are)
gCPM_imp_gene_indiv_asp_norm_sortre <- gCPM_imp_gene_indiv_asp_norm_are
head(gCPM_imp_gene_indiv_asp_norm_sortre)
rownames(gCPM_imp_gene_indiv_asp_norm_sortre) <- gCPM_imp_gene_indiv_asp_norm_sortre[,1]
gCPM_imp_gene_indiv_asp_norm_sortrebarplot <- gCPM_imp_gene_indiv_asp_norm_sortre[,-1]
gCPM_imp_gene_indiv_asp_norm_sortrebarplot <- as.matrix(t(gCPM_imp_gene_indiv_asp_norm_sortrebarplot))
head(gCPM_imp_gene_indiv_asp_norm_sortrebarplot)
dim(gCPM_imp_gene_indiv_asp_norm_sortrebarplot)
#barplot(gCPM_imp_gene_indiv_asp_norm_sortrebarplot, beside=TRUE, horiz=TRUE, col=c("red","blue","orange", "green"), xpd=FALSE,  xlim = c(0,1))


gCPM_imp_gene_indiv_asp_norm_sortrebarplot1 <- stack(gCPM_imp_gene_indiv_asp_norm_sortrebarplot)
head(gCPM_imp_gene_indiv_asp_norm_sortrebarplot1)
colnames(gCPM_imp_gene_indiv_asp_norm_sortrebarplot1) <- c("Group", "Gene", "Ratio")
gCPM_imp_gene_indiv_asp_norm_sortrebarplot1 <- data.frame(gCPM_imp_gene_indiv_asp_norm_sortrebarplot1)
head(gCPM_imp_gene_indiv_asp_norm_sortrebarplot1)
str(gCPM_imp_gene_indiv_asp_norm_sortrebarplot1)
gCPM_imp_gene_indiv_asp_norm_sortrebarplot1_half <- gCPM_imp_gene_indiv_asp_norm_sortrebarplot1
gCPM_imp_gene_indiv_asp_norm_sortrebarplot1_half["Ratio"] <- as.numeric(gCPM_imp_gene_indiv_asp_norm_sortrebarplot1_half$Ratio) - 0.5
head(gCPM_imp_gene_indiv_asp_norm_sortrebarplot1_half)
library(ggplot2)
# Basic barplot

#way1 : ggplot2
pbar <-ggplot(data=gCPM_imp_gene_indiv_asp_norm_sortrebarplot1_half, aes(x=Gene, y=Ratio, label=Group, color=Group, fill=Group)) +
  geom_bar(width = 0.5,stat="identity",  position=position_dodge(), size=0.6)
pbar + coord_flip()+
  scale_fill_manual(values=c("darkgreen","darkgreen","darkred","darkred"))+
  scale_color_manual(values=c("darkgreen","darkgreen","darkred","darkred")) + theme_bw()
ggsave("gCPM_imp_gene_indiv_asp_norm_sortrebarplot1_half.svg", width=13, height=25, units="cm", dpi=96)
ggsave("gCPM_imp_gene_indiv_asp_norm_sortrebarplot1_half.jpg", width=13, height=25, units="cm", dpi=96)
#theme(legend.position = "bottom", legend.direction = "horizontal", panel.background = element_blank(), axis.line.x = element_line(size = .6, colour = "black"))


pdot <-ggplot(data=gCPM_imp_gene_indiv_asp_norm_sortrebarplot1_half, aes(x=Gene, y=Ratio, label=Group, color="black", fill=Group)) + 
  geom_hline(yintercept = c(-0.17,0,0.17), colour = "grey", linetype=c("dashed","solid","dashed")) + 
  geom_dotplot(binaxis='y', stackdir='centerwhole', stackratio=1.5, dotsize=1.2) + theme_classic() 
pdot + coord_flip()+ ylim(c(-0.5,0.5))+
  scale_fill_manual(values=c("green","green","red","red"))+
  scale_color_manual(values=c("black","black","black","black")) 

ggsave("dotplot_gCPM_imp_gene_indiv_asp_norm_sortrebarplot1_half.svg", width=13, height=24, units="cm", dpi=96)
ggsave("dotplot_gCPM_imp_gene_indiv_asp_norm_sortrebarplot1_half.jpg", width=13, height=24, units="cm", dpi=96)



#Cluster Dotplot of 39 imprinted genes. 
#Basilia set of 37 genes  and Blcap , Adam23 added later on so 39
basilia_genes_for_clusterdotplot_as_bulk <- read.table("basilia_genes_for_clusterdotplot_as_bulk.txt", header = F)
colnames(basilia_genes_for_clusterdotplot_as_bulk) <- c("Genes")
head(gcountdataCPM_chr_gene1)
dim(gcountdataCPM_chr_gene1)
gcountdataCPM_chr_gene1_clusterdotplotgenes <- merge(basilia_genes_for_clusterdotplot_as_bulk,gcountdataCPM_chr_gene1, by ="Genes", all.x=FALSE,sort = F)
head(gcountdataCPM_chr_gene1_clusterdotplotgenes)
dim(gcountdataCPM_chr_gene1_clusterdotplotgenes)

write.csv(gcountdataCPM_chr_gene1_clusterdotplotgenes, "gcountdataCPM_chr_gene1_clusterdotplotgenes.csv", row.names = F)


#Manually done
#Open gcountdataCPM_chr_gene1_clusterdotplotgenes.csv and if 0 occurs in both allele, we replaced it with 1 (atleast one read), the other allele was still high
#save as gcountdataCPM_chr_gene1_clusterdotplotgenes_re.txt, do not copy Genes word so that header can kept true


gcountdataCPM_chr_gene1_clusterdotplotgenes_re <- read.table("gcountdataCPM_chr_gene1_clusterdotplotgenes_re.txt", header = T)
gcountdataCPM_chr_gene1_clusterdotplotgenes_re["WT_allele_ratio_1"] <- gcountdataCPM_chr_gene1_clusterdotplotgenes_re$JB1_WT_Rep1.JF1 / (gcountdataCPM_chr_gene1_clusterdotplotgenes_re$JB1_WT_Rep1.B6 + gcountdataCPM_chr_gene1_clusterdotplotgenes_re$JB1_WT_Rep1.JF1)
gcountdataCPM_chr_gene1_clusterdotplotgenes_re["WT_allele_ratio_2"] <- gcountdataCPM_chr_gene1_clusterdotplotgenes_re$JB1_WT_Rep2.JF1 / (gcountdataCPM_chr_gene1_clusterdotplotgenes_re$JB1_WT_Rep2.B6 + gcountdataCPM_chr_gene1_clusterdotplotgenes_re$JB1_WT_Rep2.JF1)
gcountdataCPM_chr_gene1_clusterdotplotgenes_re["WT_allele_ratio"] <- (gcountdataCPM_chr_gene1_clusterdotplotgenes_re$WT_allele_ratio_1 + gcountdataCPM_chr_gene1_clusterdotplotgenes_re$WT_allele_ratio_2)/2
gcountdataCPM_chr_gene1_clusterdotplotgenes_re["ZFP57KO_allele_ratio_1"] <- gcountdataCPM_chr_gene1_clusterdotplotgenes_re$JB1_ZFP57_KO_Rep1.JF1 / (gcountdataCPM_chr_gene1_clusterdotplotgenes_re$JB1_ZFP57_KO_Rep1.B6 + gcountdataCPM_chr_gene1_clusterdotplotgenes_re$JB1_ZFP57_KO_Rep1.JF1)
gcountdataCPM_chr_gene1_clusterdotplotgenes_re["ZFP57KO_allele_ratio_2"] <- gcountdataCPM_chr_gene1_clusterdotplotgenes_re$JB1_ZFP57_KO_Rep2.JF1 / (gcountdataCPM_chr_gene1_clusterdotplotgenes_re$JB1_ZFP57_KO_Rep2.B6 + gcountdataCPM_chr_gene1_clusterdotplotgenes_re$JB1_ZFP57_KO_Rep2.JF1)
gcountdataCPM_chr_gene1_clusterdotplotgenes_re["ZFP57KO_allele_ratio"] <- (gcountdataCPM_chr_gene1_clusterdotplotgenes_re$ZFP57KO_allele_ratio_1 + gcountdataCPM_chr_gene1_clusterdotplotgenes_re$ZFP57KO_allele_ratio_2)/2

head(gcountdataCPM_chr_gene1_clusterdotplotgenes_re,2)
gCPM_cluster_imp_gene_indiv_asp_norm <- gcountdataCPM_chr_gene1_clusterdotplotgenes_re

head(gCPM_cluster_imp_gene_indiv_asp_norm)
dim(gCPM_cluster_imp_gene_indiv_asp_norm)
summary(gCPM_cluster_imp_gene_indiv_asp_norm[,c(13:14,16:17)])

#Separate replicates
gCPM_cluster_imp_gene_indiv_asp_norm_are <- gCPM_cluster_imp_gene_indiv_asp_norm[,c(13:14,16:17)]
head(gCPM_cluster_imp_gene_indiv_asp_norm_are)
gCPM_cluster_imp_gene_indiv_asp_norm_sortre <- gCPM_cluster_imp_gene_indiv_asp_norm_are
head(gCPM_cluster_imp_gene_indiv_asp_norm_sortre)
gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot <- as.matrix(t(gCPM_cluster_imp_gene_indiv_asp_norm_sortre))
head(gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot)
dim(gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot)
gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1 <- stack(gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot)
head(gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1)
colnames(gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1) <- c("Group", "Gene", "Ratio")
gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1 <- data.frame(gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1)
head(gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1)
dim(gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1)
gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1["Group.1"] <- data.frame(rep(c(rep("WT_allele_ratio",2),rep("ZFP57KO_allele_ratio",2)), 39))
head(gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1)
head(gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1)
dim(gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1)
str(gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1)
gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1_half <- gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1
gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1_half["Ratio"] <- as.numeric(gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1_half$Ratio) - 0.5
head(gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1_half)
library(ggplot2)
# Dotplot for Cluster
# Change the position
pdot <-ggplot(gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1_half, aes(x=Gene, y=Ratio, fill= Group.1)) + geom_hline(yintercept = c(-0.17,0,0.17), colour = "grey", linetype=c("dashed","solid","dashed"))+
  geom_boxplot(position=position_dodge(0.8))+
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.8))
pdot + scale_fill_manual(values=c("green","red"))+ylim(c(-0.5,0.5))+
  scale_color_manual(values=c("black","black","black","black")) + theme_classic()

ggsave("dotplot_gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1_half.svg", width=90, height=10, units="cm", dpi=96)
ggsave("dotplot_gCPM_cluster_imp_gene_indiv_asp_norm_sortrebarplot1_half.jpg", width=90, height=10, units="cm", dpi=96)


#way2 : ggpubr
#library(ggpubr)
#ggbarplot(gCPM_imp_gene_indiv_asp_norm_sortrebarplot1_half, x = "Gene", y = "Ratio", fill = "Group", color = "white", palette = c("darkgreen","darkgreen","darkred","darkred"),sort.by.groups = FALSE,
#          x.text.angle = 90,          # Rotate vertically x axis texts
#          ylab = "Ratio",
#          legend.title = "Expression",
#          rotate = TRUE, position = position_dodge(),
#          ggtheme = theme_classic())
#ggdotplot(gCPM_imp_gene_indiv_asp_norm_sortrebarplot1_half, x = "Gene", y = "Ratio", fill = "Group", color = " black", palette = c("darkgreen","darkgreen","darkred","darkred"),sort.by.groups = FALSE,
#          x.text.angle = 90,          # Rotate vertically x axis texts
#          ylab = "Ratio",
#          merge = TRUE,
#          size=0.5,
#          binwidth = 0.05,
#          legend.title = "Expression",
#          rotate = TRUE,
#          repel = 2,ggsave("dotplot_gCPM_imp_gene_indiv_asp_norm_sortrebarplot1_half.jpg", width=13, height=30, units="cm", dpi=96)
#          ggtheme =  theme_classic()) + theme_cleveland() + geom_hline(yintercept = 0.0, colour = "grey") 
#ggsave("dotplot_themeclassic_gCPM_imp_gene_indiv_asp_norm_sortrebarplot1_half.svg", width=13, height=25, units="cm", dpi=96)


#theme(legend.position = "bottom", legend.direction = "horizontal", panel.background = element_blank(), axis.line.x = element_line(size = .6, colour = "black"))



#fgrep -f selected_imprinted_genes.txt /home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20_rearranged.txt -w | grep Mrpl23-ps1 -v > selected_imprinted_genes.cords.txt

#remove  Mrpl23-ps1
simprinted_gene_cords <- read.table("selected_imprinted_genes.cords.txt", header = F)
head(simprinted_gene_cords)
dim(simprinted_gene_cords)
simprinted_gene_cords <- simprinted_gene_cords[,c(1:3,5:6)]
colnames(simprinted_gene_cords) <- c("chr", "start", "end", "id", "Gene")
head(simprinted_gene_cords)
dim(simprinted_gene_cords)
head(gcountdataCPM_sortedavg_chr_gene1_ratio)
dim(gcountdataCPM_sortedavg_chr_gene1_ratio)
gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort = merge(gcountdataCPM_sortedavg_chr_gene1_ratio, simprinted_gene_cords, by="id", all.x=FALSE)
head(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort)
dim(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort)
gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_bar <- gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort[,c(5,10,11)]
head(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_bar)
library(ggplot2)
library(gplots)
rownames(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_bar)
rownames(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_bar)=gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_bar[,1]
rownames(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_bar)
colnames(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_bar)
head(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_bar)
gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot = gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_bar[,-1]
gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot = as.matrix(t(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot))
head(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot)
dim(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot)
#barplot(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot, beside=TRUE, horiz=TRUE, col=c("red","blue","orange", "green"), xpd=FALSE,  xlim = c(0,1))


gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot1 <- stack(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot)
head(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot1)
colnames(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot1) <- c("Group", "Gene", "Ratio")
gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot1 <- data.frame(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot1)
head(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot1)
str(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot1)
gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot1_half <- gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot1
gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot1_half["Ratio"] <- gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot1_half$Ratio - 0.5
head(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot1_half)
library(ggplot2)
# Basic barplot
#way1 : ggplot2
pbar <-ggplot(data=gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot1_half, aes(x=Gene, y=Ratio, label=Group, color=Group, fill=Group)) +
  geom_bar(width = 0.5,stat="identity",  position=position_dodge(), size=0.6)
pbar + coord_flip()+
  scale_fill_manual(values=c("darkgreen","darkviolet"))+
  scale_color_manual(values=c("darkgreen","darkviolet")) + theme_classic()

#Apply fisher value 0.05

gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher <- gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort
head(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher)
gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher["fisherpvalue"] <- apply(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher[,7:10],1, function(x) fisher.test(matrix(x,nrow =2))$p.value)
gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher["fisheradjpval"] <- p.adjust(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher$fisherpvalue,method="BH")
head(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher)
write.table(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher, "gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher.txt", sep="\t", quote = FALSE, append = FALSE)
gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05 <- gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher[which(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher$fisheradjpval <0.05),]
head(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05)
dim(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05)
write.table(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05, "gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05.txt", sep="\t", quote = FALSE, append = FALSE)

gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar <- gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05[,c(5,10,11)]
head(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar)
library(ggplot2)
library(gplots)
rownames(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar)
rownames(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar)=gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar[,1]
rownames(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar)
colnames(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar)
head(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar)
gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot = gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_bar[,-1]
gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot = as.matrix(t(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot))
head(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot)
dim(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot)
#barplot(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot, beside=TRUE, horiz=TRUE, col=c("red","blue","orange", "green"), xpd=FALSE,  xlim = c(0,1))


gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1 <- stack(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot)
head(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1)
colnames(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1) <- c("Group", "Gene", "Ratio")
gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1 <- data.frame(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1)
head(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1)
str(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1)
gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1_half <- gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1
gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1_half["Ratio"] <- gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1_half$Ratio - 0.5
head(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1_half)

library(ggplot2)
# Basic barplot
#way1 : ggplot2
pbar <-ggplot(data=gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1_half, aes(x=Gene, y=Ratio, label=Group, color=Group, fill=Group)) +
  geom_bar(width = 0.5,stat="identity",  position=position_dodge(), size=0.6)

pbar + coord_flip()+
  scale_fill_manual(values=c("darkgreen","darkviolet"))+
  scale_color_manual(values=c("darkgreen","darkviolet")) + theme_bw()
ggsave("gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1_half.svg", width=13, height=25, units="cm", dpi=96)
#theme(legend.position = "bottom", legend.direction = "horizontal", panel.background = element_blank(), axis.line.x = element_line(size = .6, colour = "black"))

#way2 : ggpubr
library(ggpubr)
ggbarplot(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_barplot1_half, x = "Gene", y = "Ratio", fill = "Group", color = "white", palette = c("darkgreen","darkviolet"),sort.by.groups = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Ratio",
          legend.title = "Expression",
          rotate = TRUE, position = position_dodge(),
          ggtheme = theme_classic())




ggbarplot(gcountdataCPM_sortedavg_chr_gene1_ratio_imprinted_sort_fisher_0.05_barplot1_half, x = "Gene", y = "Ratio", fill = "Group", color = "white", palette = c("darkgreen","darkviolet"),sort.by.groups = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Ratio",
          legend.title = "Expression",
          rotate = TRUE, position = position_dodge(),
          ggtheme = theme_classic())


library(WriteXLS) or use library("writexl")write_xlsx()
write_xlsx(aggregate_covdata5_per_countedbinsresep_sortbin,"recombinationmap.xlsx")
WriteXLS(gddsaSFNormcounts,ExcelFileName="AlleleSpecific_Deseq2_Normalized_Counts_all_features.xlsx", col.names=TRUE, row.names=TRUE)
WriteXLS(data.frame(z_gddsaSFNormcounts),ExcelFileName="AlleleSpecific_Deseq2_Normalized_zscores_all_features.xlsx", col.names=TRUE, row.names=TRUE)
write_xlsx(data.frame(z_gddsaSFNormcountsavg.id),"AlleleSpecific_Deseq2_Normalized_zscoresavg_all_features.xlsx")
WriteXLS(gddsaSFNormcountsavg_chr_gene1,ExcelFileName="AlleleSpecific_Deseq2_Normalized_averageCounts_all_features_gene.xlsx", col.names=TRUE, row.names=TRUE)
write_xlsx(data.frame(z_gddsaSFNormcountsavg_chr_gene1),"AlleleSpecific_Deseq2_Normalized_zscoresaverage_all_features_gene.xlsx")
WriteXLS(gddsaSFNormcounts_chr_gene1,ExcelFileName="AlleleSpecific_Deseq2_Normalized_Counts_all_features_gene.xlsx", col.names=TRUE, row.names=TRUE)
write_xlsx(data.frame(cbind(rownames(z_gddsaSFNormcountsavg.selected_genes_sort_heat1),z_gddsaSFNormcountsavg.selected_genes_sort_heat1)),"AlleleSpecific_Deseq2_zscores_selected_Imprinted_gene_heatmap.xlsx")
write_xlsx(gddsaSFNormcountsavg.selected_genes_sort_heat1,"AlleleSpecific_Deseq2_interpretation_selected_Imprinted_gene_heatmap.xlsx")
write_xlsx(gddsaSFNormcountsavg.selected_genes_sort_heat1.withratio.reorder,"AlleleSpecific_Deseq2_interpretation_selected_Imprinted_gene_updown.xlsx")
write_xlsx(z_gddsaSFNormcountsavg_imprinted_sortre,"AlleleSpecific_Deseq2_zscores_all_Imprinted_gene_heatmap.xlsx")
WriteXLS(data_imp_gene_indiv_asp_norm,ExcelFileName="AlleleSpecific_Deseq2_n_Countsandratio_imprinted_gene_verticaldotplot.xlsx", col.names=TRUE, row.names=TRUE)
WriteXLS(cluster_imp_gene_indiv_asp_norm,ExcelFileName="AlleleSpecific_Deseq2_n_Countsandratio_imprinted_gene_Clusterdotplot.xlsx", col.names=TRUE, row.names=TRUE)
WriteXLS(gcountdataCPM_sorted,ExcelFileName="AlleleSpecific_CPM_Normalized_Counts_all_features.xlsx", col.names=TRUE, row.names=TRUE)
WriteXLS(gcountdataCPM_sortedavg,ExcelFileName="AlleleSpecific_CPM_Normalized_averageCounts_all_features.xlsx", col.names=TRUE, row.names=TRUE)
WriteXLS(z_gcountdataCPM_sortedavg.id,ExcelFileName="AlleleSpecific_CPM_Normalized_zscoresavg_all_features.xlsx", col.names=TRUE, row.names=TRUE)
write_xlsx(gcountdataCPM_sortedavg_chr_gene1,"AlleleSpecific_CPM_Normalized_averageCounts_all_features_gene.xlsx")
WriteXLS(gcountdataCPM_sorted_chr_gene1,ExcelFileName="AlleleSpecific_CPM_Normalized_Counts_all_features_gene.xlsx", col.names=TRUE, row.names=TRUE)
write_xlsx(z_gcountdataCPM_sortedavg_chr_gene1,"AlleleSpecific_CPM_Normalized_zscoresaverage_all_features_gene.xlsx")
write_xlsx(data.frame(cbind(rownames(z_gcountdataCPM_sortedavg.selected_genes_sort_heat1),z_gcountdataCPM_sortedavg.selected_genes_sort_heat1)),"AlleleSpecific_CPM_zscores_selected_Imprinted_gene_heatmap.xlsx")
write_xlsx(gcountdataCPM_sortedavg.selected_genes_sort_heat1,"AlleleSpecific_CPM_interpretation_selected_Imprinted_gene_heatmap.xlsx")
write_xlsx(gcountdataCPM_sortedavg.selected_genes_sort_heat1.withratio.reorder,"AlleleSpecific_CPM_interpretation_selected_Imprinted_gene_updown.xlsx")
write_xlsx(z_gcountdataCPM_sortedavg_imprinted_sortre,"AlleleSpecific_CPM_zscores_all_Imprinted_gene_heatmap.xlsx")
WriteXLS(gCPM_imp_gene_indiv_asp_norm,ExcelFileName="AlleleSpecific_CPM_n_Countsandratio_imprinted_gene_verticaldotplot.xlsx", col.names=TRUE, row.names=TRUE)
write_xlsx(gCPM_cluster_imp_gene_indiv_asp_norm,"AlleleSpecific_CPM_n_Countsandratio_imprinted_gene_Clusterdotplot.xlsx")
write_xlsx(gddsaSFNormcountsavg_chr_gene1_fisher[,c(7,1:6)],"AlleleSpecific_gddsaSFNormcountsavg_chr_gene1_fisher.xlsx")
write_xlsx(gddsaSFNormcountsavg_chr_gene1_fisher_0.05,"AlleleSpecific_gddsaSFNormcountsavg_chr_gene1_fisher0.05.xlsx")
write_xlsx(gddsaSFNormcountsavg_chr_gene1_fisher_0.1,"AlleleSpecific_gddsaSFNormcountsavg_chr_gene1_fisher0.1.xlsx")


#Add closest distance to Deseq2 noramlized counts
grep ENS gddsaSFNormcounts_chr_gene1.txt | sort -k1,1 -k2,2n > gddsaSFNormcounts_chr_gene1.sort.txt
#While intersecting with bedtools remove duplication which happen due to overlap of Zfp57 peaks, keep only once. It is even Ok for imprinted genes even if they overlap more than one time with peaks
bedtools closest -a gddsaSFNormcounts_chr_gene1.sort.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | sort -k5,5 -u > gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10.txt
#Import this file and perform Fisher exact test and then combine
gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10 <- read.table("gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10.txt", header = F)
head(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10)
dim(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10) 
colnames(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10) <- c("chr", "start", "end", "Genes","id","JB1_WT_Rep1.B6","JB1_WT_Rep2.B6","JB1_WT_Rep1.JF1","JB1_WT_Rep2.JF1","JB1_ZFP57_KO_Rep1.B6","JB1_ZFP57_KO_Rep2.B6","JB1_ZFP57_KO_Rep1.JF1","JB1_ZFP57_KO_Rep2.JF1","chrPeak","startpeak","endPeak","peakID","qualityPeak","ClosestDistance")
head(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10)
rownames(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10) <- paste(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10$Genes,gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10$id, sep = "%")
partialgddsSFNormcounts <- gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10
#Take average
partialgddsSFNormcountsavg <- data.frame(cbind(data.frame(rownames(partialgddsSFNormcounts)),
                                               (partialgddsSFNormcounts$JB1_WT_Rep1.B6 + partialgddsSFNormcounts$JB1_WT_Rep2.B6)/2,
                                               (partialgddsSFNormcounts$JB1_WT_Rep1.JF1 + partialgddsSFNormcounts$JB1_WT_Rep2.JF1)/2, 
                                               (partialgddsSFNormcounts$JB1_ZFP57_KO_Rep1.B6 + partialgddsSFNormcounts$JB1_ZFP57_KO_Rep2.B6)/2, 
                                               (partialgddsSFNormcounts$JB1_ZFP57_KO_Rep1.JF1 + partialgddsSFNormcounts$JB1_ZFP57_KO_Rep2.JF1)/2))
head(partialgddsSFNormcountsavg)
colnames(partialgddsSFNormcountsavg) <- c("GeneID","WT.B6","WT.JF1","ZFP57KO.B6","ZFP57KO.JF1")
#Fisher exact Test
partialgddsSFNormcountsavg_fisher <- partialgddsSFNormcountsavg[,2:5]
head(partialgddsSFNormcountsavg_fisher)
partialgddsSFNormcountsavg_fisher["fisherpvalue"] <- apply(partialgddsSFNormcountsavg_fisher,1, function(x) fisher.test(matrix(x,nrow =2))$p.value)

rownames(partialgddsSFNormcountsavg_fisher) <- partialgddsSFNormcountsavg$GeneID
partialgddsSFNormcountsavg_fisher["GeneID"] <- partialgddsSFNormcountsavg$GeneID
head(partialgddsSFNormcountsavg_fisher)
#warning will be observed because we need only integers for runnung this function.
#In fisher.test(matrix(x, nrow = 2)) :'x' has been rounded to integer: Mean relative difference:

#Recombine full matrix
gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher <- cbind.data.frame(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10,partialgddsSFNormcountsavg_fisher)
head(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher,1)
gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher["fisheradjpval"] <- p.adjust(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher$fisherpvalue,method="BH")
head(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher)
write.table(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher, "gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.txt", sep="\t", quote = FALSE, append = FALSE)
#Filter p<0.05
gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher_0.05 <- gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher[which(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher$fisheradjpval <0.05),]
head(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher_0.05)
dim(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher_0.05)
write.table(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher_0.05, "gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher_0.05.txt", sep="\t", quote = FALSE, append = FALSE)

#Add strand information
awk '{if($3=="gene") print $1"\t"$4"\t"$5"\t"$3"\t"$7"\t"$10"\t"$14}' gencode.vM20.chr_patch_hapl_scaff.annotation.gtf | awk -F'"' '{print $1"\t"$2"\t"$4}' | sort -k6,6 > ens_gene_names_chrpos_dedup_M20.strand.txt
ens_gene_names_chrpos_dedup_M20.strand <- read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.strand.txt")
head(ens_gene_names_chrpos_dedup_M20.strand)
colnames(ens_gene_names_chrpos_dedup_M20.strand) <- c("chr", "start", "end","gene","strand","id","Gene")
dim(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher)
head(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher)
gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes = merge(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher, ens_gene_names_chrpos_dedup_M20.strand, by="id", all.x=TRUE)
head(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes)
gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re <- gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes[,c(27:29,31:32,1,6:24,26)]
head(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re)
dim(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re)
write_xlsx(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re, "AlleleSpecific_gddsaSFNormcounts_overlapped_ZFP57_KAP1_peaks_fisher.genes.re.xlsx")


#------------  Proportion test -------------#
#Combined p.vales -2Sumation(log(Pi))
#ChiSquare df = 2*length(p-values), so if 2 p.values 2 *2
#lower.tail =FALSE and lower.tail =TRUE, 1-p
#Prop test and tag chromsome positions
library(metaseqR)
gddsaSFNormcounts_chr_gene1 <- read.table("gddsaSFNormcounts_chr_gene1.txt", header = T)

head(gddsaSFNormcounts_chr_gene1)
dim(gddsaSFNormcounts_chr_gene1)
rownames(gddsaSFNormcounts_chr_gene1) <- paste(gddsaSFNormcounts_chr_gene1$Genes,gddsaSFNormcounts_chr_gene1$id, sep = "%")
gddsaSFNormcounts_chr_gene1_prop <- gddsaSFNormcounts_chr_gene1[,6:13]
head(gddsaSFNormcounts_chr_gene1_prop)
#Take sum of alleles
gddsaSFNormcounts_chr_gene1_prop["JB1_WT_alavg.Rep1"] <- gddsaSFNormcounts_chr_gene1_prop$JB1_WT_Rep1.B6 + gddsaSFNormcounts_chr_gene1_prop$JB1_WT_Rep1.JF1
gddsaSFNormcounts_chr_gene1_prop["JB1_WT_alavg.Rep2"] <- gddsaSFNormcounts_chr_gene1_prop$JB1_WT_Rep2.B6 + gddsaSFNormcounts_chr_gene1_prop$JB1_WT_Rep2.JF1
gddsaSFNormcounts_chr_gene1_prop["JB1_ZFP57KO_alavg.Rep1"] <- gddsaSFNormcounts_chr_gene1_prop$JB1_ZFP57_KO_Rep1.B6 + gddsaSFNormcounts_chr_gene1_prop$JB1_ZFP57_KO_Rep1.JF1
gddsaSFNormcounts_chr_gene1_prop["JB1_ZFP57KO_alavg.Rep2"] <- gddsaSFNormcounts_chr_gene1_prop$JB1_ZFP57_KO_Rep2.B6 + gddsaSFNormcounts_chr_gene1_prop$JB1_ZFP57_KO_Rep2.JF1

#Take allelic ratio
gddsaSFNormcounts_chr_gene1_prop["WT_allele_ratio_1"] <- gddsaSFNormcounts_chr_gene1_prop$JB1_WT_Rep1.JF1 / (gddsaSFNormcounts_chr_gene1_prop$JB1_WT_Rep1.B6 + gddsaSFNormcounts_chr_gene1_prop$JB1_WT_Rep1.JF1)
gddsaSFNormcounts_chr_gene1_prop["WT_allele_ratio_2"] <- gddsaSFNormcounts_chr_gene1_prop$JB1_WT_Rep2.JF1 / (gddsaSFNormcounts_chr_gene1_prop$JB1_WT_Rep2.B6 + gddsaSFNormcounts_chr_gene1_prop$JB1_WT_Rep2.JF1)
gddsaSFNormcounts_chr_gene1_prop["WT_allele_ratio"] <- (gddsaSFNormcounts_chr_gene1_prop$WT_allele_ratio_1 + gddsaSFNormcounts_chr_gene1_prop$WT_allele_ratio_2)/2
gddsaSFNormcounts_chr_gene1_prop["ZFP57KO_allele_ratio_1"] <- gddsaSFNormcounts_chr_gene1_prop$JB1_ZFP57_KO_Rep1.JF1 / (gddsaSFNormcounts_chr_gene1_prop$JB1_ZFP57_KO_Rep1.B6 + gddsaSFNormcounts_chr_gene1_prop$JB1_ZFP57_KO_Rep1.JF1)
gddsaSFNormcounts_chr_gene1_prop["ZFP57KO_allele_ratio_2"] <- gddsaSFNormcounts_chr_gene1_prop$JB1_ZFP57_KO_Rep2.JF1 / (gddsaSFNormcounts_chr_gene1_prop$JB1_ZFP57_KO_Rep2.B6 + gddsaSFNormcounts_chr_gene1_prop$JB1_ZFP57_KO_Rep2.JF1)
gddsaSFNormcounts_chr_gene1_prop["ZFP57KO_allele_ratio"] <- (gddsaSFNormcounts_chr_gene1_prop$ZFP57KO_allele_ratio_1 + gddsaSFNormcounts_chr_gene1_prop$ZFP57KO_allele_ratio_2)/2

dim(gddsaSFNormcounts_chr_gene1_prop) #15673    18
gddsaSFNormcounts_chr_gene1_propid <- gddsaSFNormcounts_chr_gene1_prop
gddsaSFNormcounts_chr_gene1_propid["id"] <- rownames(gddsaSFNormcounts_chr_gene1_prop)
head(gddsaSFNormcounts_chr_gene1_propid)
writexl::write_xlsx(gddsaSFNormcounts_chr_gene1_propid, "gddsaSFNormcounts_chr_gene1_prop.xlsx")
#Filter B6+JF1 > 10 in WT, if both the alleles are 10 in either of the replicates in WT
gddsaSFNormcounts_chr_gene1_prop.filt <- gddsaSFNormcounts_chr_gene1_prop[which(gddsaSFNormcounts_chr_gene1_prop$JB1_WT_alavg.Rep1 > 10 &
                                                                                  gddsaSFNormcounts_chr_gene1_prop$JB1_WT_alavg.Rep2 > 10),]


dim(gddsaSFNormcounts_chr_gene1_prop.filt)#11371    18
head(gddsaSFNormcounts_chr_gene1_prop.filt,1)
summary(gddsaSFNormcounts_chr_gene1_prop.filt)
#WT
JB1_WT_alavg.Rep1.prop1 <- Map(prop.test,x =gddsaSFNormcounts_chr_gene1_prop.filt$JB1_WT_Rep1.B6, n= gddsaSFNormcounts_chr_gene1_prop.filt$JB1_WT_alavg.Rep1, p=0.5)

gddsaSFNormcounts_chr_gene1_prop.filt["JB1_WT_alavg.Rep1.pvalue"] <- data.frame(capture.output(for (i in seq_along(JB1_WT_alavg.Rep1.prop1)){
  cat(JB1_WT_alavg.Rep1.prop1[[i]]$p.value, "\n")
}))

JB1_WT_alavg.Rep2.prop1 <- Map(prop.test,x =gddsaSFNormcounts_chr_gene1_prop.filt$JB1_WT_Rep2.B6, n= gddsaSFNormcounts_chr_gene1_prop.filt$JB1_WT_alavg.Rep2, p=0.5)

gddsaSFNormcounts_chr_gene1_prop.filt["JB1_WT_alavg.Rep2.pvalue"] <- data.frame(capture.output(for (i in seq_along(JB1_WT_alavg.Rep2.prop1)){
  cat(JB1_WT_alavg.Rep2.prop1[[i]]$p.value, "\n")
}))

#Combine p-values calculated from prop.test for rep1 and rep2
gddsaSFNormcounts_chr_gene1_prop.filt["JB1_WT_alavg.comb.testpvalue"] <- -2 *(log(as.numeric(as.character(gddsaSFNormcounts_chr_gene1_prop.filt$JB1_WT_alavg.Rep1.pvalue))) + log(as.numeric(as.character(gddsaSFNormcounts_chr_gene1_prop.filt$JB1_WT_alavg.Rep2.pvalue))))
#To deal with -2 * (log(0)+log(0)), I used fisher.method function which has option to zero.sub to value, I set it 2.2e-16, this will replace Inf with 2.2e-16, 
#eg. -2*(log(2.2e-16)+log(2.2e-16)) = 144.2116
class(gddsaSFNormcounts_chr_gene1_prop.filt)
gddsaSFNormcounts_chr_gene1_prop.filt <- data.frame(gddsaSFNormcounts_chr_gene1_prop.filt)
gddsaSFNormcounts_chr_gene1_prop.filt["JB1_WT_alavg.comb.fisherpvalue"] <- fisher.method(data.frame(as.numeric(as.character(gddsaSFNormcounts_chr_gene1_prop.filt$JB1_WT_alavg.Rep1.pvalue)), as.numeric(as.character(gddsaSFNormcounts_chr_gene1_prop.filt$JB1_WT_alavg.Rep2.pvalue))), method = c("fisher"), p.corr ="none", zero.sub = 2.2e-16, na.rm = FALSE, mc.cores=NULL)
#Warning will pop up Warning message:
#In `[<-.data.frame`(`*tmp*`, "JB1_WT_alavg.comb.fisherpvalue", value = list( :provided 4 variables to replace 1 variables
#Export the sheet and manually checked the column JB1_WT_alavg.comb.testpvalue and JB1_WT_alavg.comb.fisherpvalue should match except Inf, so manual and fisher.methhod is ok
gddsaSFNormcounts_chr_gene1_prop.filt_wt <- gddsaSFNormcounts_chr_gene1_prop.filt
gddsaSFNormcounts_chr_gene1_prop.filt_wt["ID"] <- rownames(gddsaSFNormcounts_chr_gene1_prop.filt_wt)
writexl::write_xlsx(gddsaSFNormcounts_chr_gene1_prop.filt_wt, "gddsaSFNormcounts_chr_gene1_prop.filt_check1forWT.xlsx")

head(gddsaSFNormcounts_chr_gene1_prop.filt)
tail(gddsaSFNormcounts_chr_gene1_prop.filt)
#Chisquare test for p-values#pchisq(q, df, ncp = 0, lower.tail = TRUE, log.p = FALSE)#p vector of probabilities, n	 number of observations, degrees of freedom (non-negative, but can be non-integer)
gddsaSFNormcounts_chr_gene1_prop.filt["JB1_WT_alavg.comb.pvalue.pchisq"] <- pchisq(gddsaSFNormcounts_chr_gene1_prop.filt$JB1_WT_alavg.comb.fisherpvalue,4, lower.tail=FALSE)
gddsaSFNormcounts_chr_gene1_prop.filt["JB1_WT_alavg.comb.fdr"] <- p.adjust(gddsaSFNormcounts_chr_gene1_prop.filt$JB1_WT_alavg.comb.pvalue.pchisq,method="BH")
head(gddsaSFNormcounts_chr_gene1_prop.filt)



#Introduce gene name
gddsaSFNormcounts_chr_gene1_prop.filt["id"] <- data.frame(rownames(gddsaSFNormcounts_chr_gene1_prop.filt))
head(gddsaSFNormcounts_chr_gene1_prop.filt)
library(splitstackshape)
gddsaSFNormcounts_chr_gene1_prop.filt.sep <- cSplit(gddsaSFNormcounts_chr_gene1_prop.filt, "id", "%")
dim(gddsaSFNormcounts_chr_gene1_prop.filt.sep)
head(gddsaSFNormcounts_chr_gene1_prop.filt.sep)
gddsaSFNormcounts_chr_gene1_prop.filt_gene = data.frame(gddsaSFNormcounts_chr_gene1_prop.filt.sep[,c(25,26,1:24)])
head(gddsaSFNormcounts_chr_gene1_prop.filt_gene)
write.table(gddsaSFNormcounts_chr_gene1_prop.filt_gene, "gddsaSFNormcounts_chr_gene1_prop.filt_gene.txt", sep = "\t", quote = F, row.names = F, append = F)
writexl::write_xlsx(gddsaSFNormcounts_chr_gene1_prop.filt_gene,"AlleleSpecific_gddsaSFNormcounts_chr_gene1_proptest_all.xlsx")


#Assign chr to gddsaSFNormcounts_chr_gene1_prop.filt_gene
chr.position  =  read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt",header=FALSE)
head(chr.position)
colnames(chr.position) <- c("id_2", "Genes", "chr", "start", "end")
head(chr.position)
dim(chr.position)
head(gddsaSFNormcounts_chr_gene1_prop.filt_gene)
dim(gddsaSFNormcounts_chr_gene1_prop.filt_gene)
gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr = merge(gddsaSFNormcounts_chr_gene1_prop.filt_gene, chr.position, by="id_2", all.x=FALSE)
head(gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr)
dim(gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr)
gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr <- gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr[,c(27:30,1:26)]
write.table(gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr, "gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr.txt", sep = "\t", quote = F, row.names = F, append = F)


#Monoalleically expressed in Wildtype
gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig <- gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr[which(gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr$JB1_WT_alavg.comb.fdr<0.05),]
head(gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig)
dim(gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig)
write.table(gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig, "gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig.txt", sep = "\t", quote = F, row.names = F, append = F)
writexl::write_xlsx(gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig,"AlleleSpecific_gddsaSFNormcounts_chr_gene1_proptest_WT_fdrsig.xlsx")
gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig0.1 <- gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr[which(gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr$JB1_WT_alavg.comb.fdr<0.1),]
head(gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig0.1)
dim(gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig0.1)
write.table(gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig0.1, "gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig0.1.txt", sep = "\t", quote = F, row.names = F, append = F)
writexl::write_xlsx(gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig0.1,"AlleleSpecific_gddsaSFNormcounts_chr_gene1_proptest_WT_fdrsig0.1.xlsx")


#Assign Imprinted gene names
fgrep -f imprinted_gene_name.txt gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig.txt -w > imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig.txt 
#remove chrX, (chrY, chrM) already not present in  imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig.txt
grep chrX imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig.txt -v > imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX.txt
imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX <- read.table("imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX.txt", header = F, stringsAsFactors = F)
head(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX)
colnames(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX) <- colnames(gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig)
head(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX)
dim(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX)
writexl::write_xlsx(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX,"imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX.xlsx")

#Add id as it is required for fisher merge
imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX["id"] <- imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX$id_2
dim(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX)
head(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX)


gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher <- read.table("gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.txt", header = T)
#Add strand information
awk '{if($3=="gene") print $1"\t"$4"\t"$5"\t"$3"\t"$7"\t"$10"\t"$14}' gencode.vM20.chr_patch_hapl_scaff.annotation.gtf | awk -F'"' '{print $1"\t"$2"\t"$4}' | sort -k6,6 > ens_gene_names_chrpos_dedup_M20.strand.txt
ens_gene_names_chrpos_dedup_M20.strand <- read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.strand.txt")
head(ens_gene_names_chrpos_dedup_M20.strand)
colnames(ens_gene_names_chrpos_dedup_M20.strand) <- c("chr", "start", "end","gene","strand","id","Gene")
dim(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher)
head(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher)
gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes = merge(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher, ens_gene_names_chrpos_dedup_M20.strand, by="id", all.x=TRUE)
head(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes)
gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re <- gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes[,c(27:29,31:32,1,6:24,26)]
head(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re)
dim(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re)



imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher <- merge(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re, imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX, by="id", all.x=F)
head(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher,1)
dim(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher)
writexl::write_xlsx(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher,"imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher.xlsx")

#Note in this data.frame fisherpvalue fisheradjpval are allelic disbalance 2x2 matrix fisher test while *.comb.fisherpvalue are from prop test p-value combine, They are different
#Now all imprinted genes below are statisitcal significant from Prop test in WT
#Filter by allelic ratio
#rep1
imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af <- imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher[which(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher$WT_allele_ratio_1 <= 0.33 |
                                                                                                                                                                                                                                      imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher$WT_allele_ratio_1 >= 0.67),]

dim(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af)
#rep2
imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af2 <- imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af[which(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af$WT_allele_ratio_2 <= 0.33 |
                                                                                                                                                                                                                                          imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af$WT_allele_ratio_2 >= 0.67),]
dim(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af2)
#Remove very low expresed genes
imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af3 <- imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af2[which(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af2$Gene != "B830012L14Rik"),]
#I also Remove Th as it is on boundary line in fdr and very low expressed
imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af3 <- imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af3[which(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af3$Gene != "Th"),]

dim(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af3)
#See imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher.xlsx for manual filtering 32 were left
head(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af3)
writexl::write_xlsx(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af3,
           "imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af3.xlsx")




#Do it for all prop genes
#Add id as it is required for fisher merge
gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr["id"] <- gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr$id_2
dim(gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr)
head(gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr)

head(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re)
gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher <- merge(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re, gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr, by="id", all.x=TRUE)
head(gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher,1)
dim(gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher)
writexl::write_xlsx(gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher,"gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher.xlsx")
#Gene with NA in prop data column are not covered in fisher test, as checked manually




#Re-Dotplot, only significant prop test and sufficiently expressed,  imprinted genes
data_for_merge_dotplot_NPCs <- imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af3[,c(6,47,50,56)]
head(data_for_merge_dotplot_NPCs)
colnames(data_for_merge_dotplot_NPCs) <- c("NPCs_Genes", "NPCs_WT_allele_ratio", "NPCs_ZFP57KO_allele_ratio", "NPCs_JB1_WT_alavg.comb.fdr")
data_for_merge_dotplot_NPCs["Genes"] <- data_for_merge_dotplot_NPCs$NPCs_Genes
write.table(data_for_merge_dotplot_NPCs$NPCs_Genes, "data_for_merge_dotplot_NPCs_Genes_sigsel.txt", col.names = F,row.names = F, quote = F, append = F)

data_for_merge_dotplot_ESCs <- read.table("/media/ankitv/Archivio1/2021/alsp_mm10_T0/featurecounts/Igf2/data_for_merge_dotplot_ESCs.txt", header = T)
head(data_for_merge_dotplot_ESCs)
data_for_merge_dotplot_ESCs_NPCs = merge(data_for_merge_dotplot_ESCs,data_for_merge_dotplot_NPCs, by.x="Genes",by.y="Genes", all.x=F, all.y=TRUE)
head(data_for_merge_dotplot_ESCs_NPCs)
dim(data_for_merge_dotplot_ESCs_NPCs)
writexl::write_xlsx(data_for_merge_dotplot_ESCs_NPCs,"data_for_merge_dotplot_ESCs_NPCs_withFDR.xlsx")
rownames(data_for_merge_dotplot_ESCs_NPCs) <- data_for_merge_dotplot_ESCs_NPCs$Genes
data_for_merge_dotplot_ESCs_NPCs <- data_for_merge_dotplot_ESCs_NPCs[,-1]
head(data_for_merge_dotplot_ESCs_NPCs)

data_for_merge_dotplot_ESCs_NPCs <- data_for_merge_dotplot_ESCs_NPCs[,c(2,3,6,7)]
head(data_for_merge_dotplot_ESCs_NPCs)
#Sort Monoalleic in NPC_WT
data_for_merge_dotplot_ESCs_NPCs <- data_for_merge_dotplot_ESCs_NPCs[order(data_for_merge_dotplot_ESCs_NPCs$NPCs_WT_allele_ratio),]
data_for_merge_dotplot_ESCs_NPCs <- as.matrix(t(data_for_merge_dotplot_ESCs_NPCs))
head(data_for_merge_dotplot_ESCs_NPCs)
dim(data_for_merge_dotplot_ESCs_NPCs)

#barplot(data_for_merge_dotplot_ESCs_NPCs, beside=TRUE, horiz=TRUE, col=c("red","blue","orange", "green"), xpd=FALSE,  xlim = c(0,1))
data_for_merge_dotplot_ESCs_NPCst <- stack(data_for_merge_dotplot_ESCs_NPCs)
head(data_for_merge_dotplot_ESCs_NPCst)

colnames(data_for_merge_dotplot_ESCs_NPCst) <- c("Group", "Gene", "Ratio")
data_for_merge_dotplot_ESCs_NPCst <- data.frame(data_for_merge_dotplot_ESCs_NPCst)
head(data_for_merge_dotplot_ESCs_NPCst)
str(data_for_merge_dotplot_ESCs_NPCst)
data_for_merge_dotplot_ESCs_NPCst_half <- data_for_merge_dotplot_ESCs_NPCst
data_for_merge_dotplot_ESCs_NPCst_half["Ratio"] <- as.numeric(data_for_merge_dotplot_ESCs_NPCst_half$Ratio) - 0.5
head(data_for_merge_dotplot_ESCs_NPCst_half)
write.table(data_for_merge_dotplot_ESCs_NPCst_half, "data_for_merge_dotplot_ESCs_NPCst_half.txt", sep = "\t", quote = F, append = F)

psigescnpc <-ggplot(data=data_for_merge_dotplot_ESCs_NPCst_half, aes(x=Gene, y=Ratio, label=Group, color="black", fill=Group)) + 
  geom_hline(yintercept = c(-0.17,0,0.17), colour = "grey", linetype=c("dashed","solid","dashed")) + 
  geom_point(aes(shape=Group, color=Group, fill=Group, size=Group), position = position_dodge(width = 0.1)) + theme_classic()  

psigescnpc + coord_flip()+ ylim(c(-0.5,0.5))+
  scale_fill_manual(values=c("#3CB371","#B22222","green","red"))+
  scale_color_manual(values=c("black","black","black","black")) +
  scale_shape_manual(values=c(22,22,24,24))+
  scale_size_manual(values=c(2,2,2,2)) 

ggsave("psigescnpcplot_data_imp_gene_indiv_asp_norm_resig_sortrebarplot1_half.svg", width=14, height=20, units="cm", dpi=96)
ggsave("psigescnpcplot_data_imp_gene_indiv_asp_norm_resig_sortrebarplot1_half.jpg", width=14, height=20, units="cm", dpi=96)

#Tablewithchr
data_for_merge_dotplot_NPCschr <- imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af3[,c(2:4,1,5,6,47,50,56)]
head(data_for_merge_dotplot_NPCschr)
colnames(data_for_merge_dotplot_NPCschr) <- c("chr", "start", "end","ensid","strand","NPCs_Genes", "NPCs_WT_allele_ratio", "NPCs_ZFP57KO_allele_ratio", "NPCs_JB1_WT_alavg.comb.fdr")
data_for_merge_dotplot_NPCschr["Genes"] <- data_for_merge_dotplot_NPCschr$NPCs_Genes

data_for_merge_dotplot_ESCs_NPCstab = merge(data_for_merge_dotplot_ESCs,data_for_merge_dotplot_NPCschr, by.x="Genes",by.y="Genes", all.x=F, all.y=TRUE)
head(data_for_merge_dotplot_ESCs_NPCstab)
dim(data_for_merge_dotplot_ESCs_NPCstab)
data_for_merge_dotplot_ESCs_NPCstab <- data_for_merge_dotplot_ESCs_NPCstab[,c(6:8,9,10,1,2:5,11:14)]
head(data_for_merge_dotplot_ESCs_NPCstab,2)
data_for_merge_dotplot_ESCs_NPCstab <- data_for_merge_dotplot_ESCs_NPCstab[order(data_for_merge_dotplot_ESCs_NPCstab$chr, data_for_merge_dotplot_ESCs_NPCstab$start),]
writexl::write_xlsx(data_for_merge_dotplot_ESCs_NPCstab,"Allele_specific_data_for_merge_dotplot_ESCs_NPCstab_withFDR.xlsx")

#NPC significant and ESC sufficiently expressed
selectednpcsig_dotplot_ESCs <- read.table("/media/ankitv/Archivio1/2021/alsp_mm10_T0/featurecounts/Igf2/selectednpcsig_dotplot_ESCs.txt", header=T)
#Get only NPC sig genes
selectednpcsig_dotplot_ESCs_NPCs = merge(selectednpcsig_dotplot_ESCs,data_for_merge_dotplot_NPCs, by.x="Genes",by.y="Genes", all.x=F, all.y=TRUE)
head(selectednpcsig_dotplot_ESCs_NPCs)
dim(selectednpcsig_dotplot_ESCs_NPCs)
writexl::write_xlsx(selectednpcsig_dotplot_ESCs_NPCs,"selectednpcsig_dotplot_ESCs_NPCs_withFDR.xlsx")
rownames(selectednpcsig_dotplot_ESCs_NPCs) <- selectednpcsig_dotplot_ESCs_NPCs$Genes
selectednpcsig_dotplot_ESCs_NPCs <- selectednpcsig_dotplot_ESCs_NPCs[,-1]
head(selectednpcsig_dotplot_ESCs_NPCs)

selectednpcsig_dotplot_ESCs_NPCs <- selectednpcsig_dotplot_ESCs_NPCs[,c(2,3,5,6)]
head(selectednpcsig_dotplot_ESCs_NPCs)
#Sort Monoalleic in NPC_WT
selectednpcsig_dotplot_ESCs_NPCs <- selectednpcsig_dotplot_ESCs_NPCs[order(selectednpcsig_dotplot_ESCs_NPCs$NPCs_WT_allele_ratio),]
selectednpcsig_dotplot_ESCs_NPCs <- as.matrix(t(selectednpcsig_dotplot_ESCs_NPCs))
head(selectednpcsig_dotplot_ESCs_NPCs)
dim(selectednpcsig_dotplot_ESCs_NPCs)

#barplot(selectednpcsig_dotplot_ESCs_NPCs, beside=TRUE, horiz=TRUE, col=c("red","blue","orange", "green"), xpd=FALSE,  xlim = c(0,1))
selectednpcsig_dotplot_ESCs_NPCst <- stack(selectednpcsig_dotplot_ESCs_NPCs)
head(selectednpcsig_dotplot_ESCs_NPCst)

colnames(selectednpcsig_dotplot_ESCs_NPCst) <- c("Group", "Gene", "Ratio")
selectednpcsig_dotplot_ESCs_NPCst <- data.frame(selectednpcsig_dotplot_ESCs_NPCst)
head(selectednpcsig_dotplot_ESCs_NPCst)
str(selectednpcsig_dotplot_ESCs_NPCst)
selectednpcsig_dotplot_ESCs_NPCst_half <- selectednpcsig_dotplot_ESCs_NPCst
selectednpcsig_dotplot_ESCs_NPCst_half["Ratio"] <- as.numeric(selectednpcsig_dotplot_ESCs_NPCst_half$Ratio) - 0.5
head(selectednpcsig_dotplot_ESCs_NPCst_half)
write.table(selectednpcsig_dotplot_ESCs_NPCst_half, "selectednpcsig_dotplot_ESCs_NPCst_half.txt", sep = "\t", quote = F, append = F)

pescsignpc <-ggplot(data=selectednpcsig_dotplot_ESCs_NPCst_half, aes(x=Gene, y=Ratio, label=Group, color="black", fill=Group)) + 
  geom_hline(yintercept = c(-0.17,0,0.17), colour = "grey", linetype=c("dashed","solid","dashed")) + 
  geom_point(aes(shape=Group, color=Group, fill=NA, size=Group), stroke = 2, position = position_dodge(width = 0.0)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA),  plot.background = element_rect(fill = "transparent",colour = NA)) 


pescsignpc + coord_flip()+ ylim(c(-0.5,0.5))+
  scale_fill_manual(values=c("#00cc00","#cd0000","green","red"))+
  scale_color_manual(values=c("#00cc00","#cd0000","green","red")) +
  scale_shape_manual(values=c(21,21,24,24))+
  scale_size_manual(values=c(4,4,3.5,3.5)) 
#Warning message:Removed 26 rows containing missing values (geom_point).# 13 NA for Wt and 13 NA for ZFP57KO in ESCs.
ggsave("pescsignpcplot_data_imp_gene_indiv_asp_norm_resig_sortrebarplot1_half.png", width=22, height=48, units="cm", dpi=96,bg = "transparent")
ggsave("pescsignpcplot_data_imp_gene_indiv_asp_norm_resig_sortrebarplot1_half.svg", width=22, height=48, units="cm", dpi=96,bg = "transparent")
ggsave("pescsignpcplot_data_imp_gene_indiv_asp_norm_resig_sortrebarplot1_half.jpg", width=22, height=48, units="cm", dpi=96,bg = "transparent")

#Tablewithchr
data_for_merge_dotplot_NPCsrchr <- imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af3[,c(2:4,1,5,6,45:50,56)]
head(data_for_merge_dotplot_NPCsrchr)
colnames(data_for_merge_dotplot_NPCsrchr) <- c("chr", "start", "end","ensid","strand","NPCs_Genes", "NPCs_WT_allele_ratio_1","NPCs_WT_allele_ratio_2","NPCs_WT_allele_ratio", "NPCs_ZFP57KO_allele_ratio_1", "NPCs_ZFP57KO_allele_ratio_2", "NPCs_ZFP57KO_allele_ratio", "NPCs_JB1_WT_alavg.comb.fdr")
data_for_merge_dotplot_NPCsrchr["Genes"] <- data_for_merge_dotplot_NPCsrchr$NPCs_Genes
dim(data_for_merge_dotplot_NPCsrchr)

selectednpcsig_dotplot_ESCs_NPCstab = merge(selectednpcsig_dotplot_ESCs,data_for_merge_dotplot_NPCsrchr, by.x="Genes",by.y="Genes", all.x=F, all.y=TRUE)
head(selectednpcsig_dotplot_ESCs_NPCstab)
dim(selectednpcsig_dotplot_ESCs_NPCstab)
selectednpcsig_dotplot_ESCs_NPCstab <- selectednpcsig_dotplot_ESCs_NPCstab[,c(5:7,8,9,1,2:4,10:17)]
head(selectednpcsig_dotplot_ESCs_NPCstab,2)
selectednpcsig_dotplot_ESCs_NPCstab <- selectednpcsig_dotplot_ESCs_NPCstab[order(selectednpcsig_dotplot_ESCs_NPCstab$chr, selectednpcsig_dotplot_ESCs_NPCstab$start),]
writexl::write_xlsx(selectednpcsig_dotplot_ESCs_NPCstab,"Allele_specific_selectednpcsig_dotplot_ESCs_NPCstab_withFDR.xlsx")


#Prepare Alsp Table all remianing genes after filtered and prop.test 
gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher_shared <- merge(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re, gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr, by="id", all.y=TRUE, all.x=FALSE)
head(gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher_shared,1)
dim(gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher_shared)
Allele_Specific_Supplementary_Table_4  <- gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher_shared
Allele_Specific_Supplementary_Table_4 <- Allele_Specific_Supplementary_Table_4[,c(2:6,1,7:24,45:50,56)]
head(Allele_Specific_Supplementary_Table_4)
colnames(Allele_Specific_Supplementary_Table_4) <- c("chr","start","end","strand","Gene","id","JB1_WT_Rep1.B6","JB1_WT_Rep2.B6","JB1_WT_Rep1.JF1","JB1_WT_Rep2.JF1","JB1_ZFP57_KO_Rep1.B6","JB1_ZFP57_KO_Rep2.B6","JB1_ZFP57_KO_Rep1.JF1","JB1_ZFP57_KO_Rep2.JF1","chrPeak","startpeak","endPeak","peakID","qualityPeak","ClosestDistance","WT.B6","WT.JF1","ZFP57KO.B6","ZFP57KO.JF1","WT_allele_ratio_1","WT_allele_ratio_2","WT_allele_ratio","ZFP57KO_allele_ratio_1","ZFP57KO_allele_ratio_2","ZFP57KO_allele_ratio","JB1_WT_alavg.comb.fdr")
write.table(Allele_Specific_Supplementary_Table_4, "Allele_Specific_Supplementary_Table_4.txt", sep="\t", row.names = F, quote = F, append = F)
writexl::write_xlsx(Allele_Specific_Supplementary_Table_4,"Allele_Specific_Supplementary_Table_4.xlsx")

#Prepare Alsp Table all genes in ASpecfic genes 15673  genes data and prop.test 
head(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re)
colnames(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re) <- c("chr","start","end","strand","Gene","id","JB1_WT_Rep1.B6","JB1_WT_Rep2.B6","JB1_WT_Rep1.JF1","JB1_WT_Rep2.JF1","JB1_ZFP57_KO_Rep1.B6","JB1_ZFP57_KO_Rep2.B6","JB1_ZFP57_KO_Rep1.JF1","JB1_ZFP57_KO_Rep2.JF1","chrPeak","startpeak","endPeak","peakID","qualityPeak","ClosestDistance","WT.B6","WT.JF1","ZFP57KO.B6","ZFP57KO.JF1","fisherpvalue", "fisheradjpval")

#Take allelic ratio
gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio <- gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re
gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio["WT_allele_ratio_1"] <- gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio$JB1_WT_Rep1.JF1 / (gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio$JB1_WT_Rep1.B6 + gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio$JB1_WT_Rep1.JF1)
gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio["WT_allele_ratio_2"] <- gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio$JB1_WT_Rep2.JF1 / (gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio$JB1_WT_Rep2.B6 + gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio$JB1_WT_Rep2.JF1)
gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio["WT_allele_ratio"] <- (gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio$WT_allele_ratio_1 + gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio$WT_allele_ratio_2)/2
gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio["ZFP57KO_allele_ratio_1"] <- gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio$JB1_ZFP57_KO_Rep1.JF1 / (gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio$JB1_ZFP57_KO_Rep1.B6 + gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio$JB1_ZFP57_KO_Rep1.JF1)
gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio["ZFP57KO_allele_ratio_2"] <- gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio$JB1_ZFP57_KO_Rep2.JF1 / (gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio$JB1_ZFP57_KO_Rep2.B6 + gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio$JB1_ZFP57_KO_Rep2.JF1)
gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio["ZFP57KO_allele_ratio"] <- (gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio$ZFP57KO_allele_ratio_1 + gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio$ZFP57KO_allele_ratio_2)/2
head(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio)
gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher_shared2 <- merge(gddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_withratio,Allele_Specific_Supplementary_Table_4 , by="id", all.y=TRUE, all.x=TRUE)
head(gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher_shared2,2)
tail(gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher_shared2)
dim(gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher_shared2)
write.table(gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher_shared2, "gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher_shared2.txt", sep="\t", row.names = F, quote = F, append = F)

Allele_Specific_ALL_Supplementary_Table_4  <- gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher_shared2[order(gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher_shared2$chr.x, gddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher_shared2$start.x),]
#take chr gene start end starnd info and ratio  for all genes
colnames(Allele_Specific_ALL_Supplementary_Table_4[,c(2:6,1,7:20,27:32,62)])
Allele_Specific_ALL_Supplementary_Table_4 <- Allele_Specific_ALL_Supplementary_Table_4[,c(2:6,1,7:20,27:32,62)]
head(Allele_Specific_ALL_Supplementary_Table_4)
tail(Allele_Specific_ALL_Supplementary_Table_4)
write.table(Allele_Specific_ALL_Supplementary_Table_4, "Allele_Specific_ALL_Supplementary_Table_4.txt", sep="\t", row.names = F, quote = F, append = F)

#Coloring
fgrep -f imprinted_gene_name.txt  Allele_Specific_ALL_Supplementary_Table_4.txt -w > Allele_Specific_ALL_Supplementary_Table_4_imprinted_genes.txt
fgrep -f imprinted_gene_name.txt  Allele_Specific_ALL_Supplementary_Table_4.txt -w -v > Allele_Specific_ALL_Supplementary_Table_4_Nonimprinted_genes.txt
#Open one by one and paste in excel sheet
#First paste imprinted color red all genes
#Now paste non-imprinted and take header and put it in first line
#Benefit of doing this it will also take care of NA and mention blank as NA
#Color by alleles, sort by chr and start---Manually checked the sorting if done proterly
#Name columns as required, 
#Also in excel sheet manually filtered NA from adj p-value, leftover is 11217 which is equal to dataframe used for prop test
#data in excel test_roughNA with NA shows either B6+JF1 <10 in WT or B6+JF1 <0 in ZFP57KO
#save as Allele_Specific_ALL_Supplementary_Table_4.xlsx. 
#Divide the sheet in two parts NA and without NA

########################  End of Analysis ##################################à

