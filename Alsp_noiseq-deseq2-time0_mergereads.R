########################################### Featurecounts ####################################################

#igf2, SE gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschrlike-.gtf can be used as it has no readlength isssue
##/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-GRCm38-M20-overlapped/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschrlike-.gtf -T 12 -o Alsp_JB1_time0_star-featureCounts_GRCm38.mm10_Igf2.txt JB1_WT_Aligned.sortedByReadname.out.genome.sort.B6.bam  JB1_WT_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Aligned.sortedByReadname.out.genome.sort.JF1.bam

################################################################# DESeq2 ######################################################################
setwd("/media/ankitv/Archivio1/2021/alsp_mm10_T0/featurecounts/Igf2")
countdata <- read.table("Alsp_JB1_time0_star-featureCounts_GRCm38.mm10_Igf2.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]
head(countdata)
# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))
colnames(countdata) <-  c("JB1_WT_B6","JB1_WT_JF1", "JB1_ZFP57KO_B6", "JB1_ZFP57KO_JF1")
head(countdata)

###################################################### DESeq2 #######################################################
library(DESeq2)
countData <- as.matrix(countdata)
coldata <- read.table("coldata.txt" , header = TRUE, stringsAsFactors = FALSE)
rownames(coldata)
coldata[,1]
rownames(coldata)=coldata[,1]
rownames(coldata)
colnames(coldata)
coldata = coldata[,-1]
coldata <- data.frame(coldata)
head(coldata)
coldata <- coldata[,c("condition","replicate","allele")]
coldata$condition <- factor(coldata$condition)
rownames(coldata) <- colnames(countData)
all(rownames(coldata) == colnames(countData)) #should print TRUE
dds <- DESeqDataSetFromMatrix(countData =countData, colData = coldata, design = ~ condition)
dds
#featureData <- data.frame(gene=rownames(countData))
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#View filtered count matrix: View(counts(dds))
#Apply first method Size factor based of DESeq2
ddsSF <- dds
ddsSF <- estimateSizeFactors(ddsSF)
sizeFactors(ddsSF)
#Note normalization =FALSE will not divide counts by normalization factors
ddsSFcounts <- counts(ddsSF, normalized=FALSE)
head(ddsSFcounts)
boxplot(ddsSFcounts, ylim=c(0,200))
#Normalization is the part of DESeq command: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
ddsSFNorm <- dds
#Bulk : 0.8250728       1.2120143 of Combinefersmithlab
#sizeFactors(ddsSFNorm) <- c(0.8250728, 0.8250728, 1.2120143, 1.2120143)
#Bulk : c(0.8248221, 1.2123827) of internal fastq
sizeFactors(ddsSFNorm) <- c(0.8248221, 0.8248221, 1.2123827, 1.2123827)
sizeFactors(ddsSFNorm)

#Normalize by DESeq2 size factor, set size factoe using Bdds 
#Note normalization =TRUE divide counts by the user set normalization factors
#Normalize allele specfic counts export normalized counts
ddsSFNormcounts <- counts(ddsSFNorm, normalized=TRUE)
head(ddsSFNormcounts)
write.table(ddsSFNormcounts, "ddsSFNormcounts.txt", sep="\t", quote=F, col.names=NA)
#DEG
#With no replicates deseq will not work. Since I only needed normalized count matrix I assigned gene name and took the matrix out
#Chromosome positions
chr.pos =  read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt",header=FALSE)
head(chr.pos)
colnames(chr.pos) <- c("ensid", "Genes", "chr", "start", "end")
head(chr.pos)
head(ddsSFNormcounts)
ddsSFNormcounts <- data.frame(ddsSFNormcounts)
ddsSFNormcounts["ensid"] <- rownames(ddsSFNormcounts)
head(ddsSFNormcounts)
ddsSFNormcounts_chr.pos = merge(ddsSFNormcounts, chr.pos, by="ensid", all.x=TRUE)
head(ddsSFNormcounts_chr.pos)
dim(ddsSFNormcounts_chr.pos)
ddsSFNormcounts_chr.pos <- ddsSFNormcounts_chr.pos[,c(7:9,6,1:5)]
write.table(ddsSFNormcounts_chr.pos, "ddsSFNormcounts_chr.pos_gene.txt", sep="\t", quote = FALSE, append = FALSE)
boxplot(ddsSFNormcounts[,1:4], ylim=c(0,500))

ddsaSFNormcounts_chr_gene1 <- ddsSFNormcounts_chr.pos
head(ddsaSFNormcounts_chr_gene1)
write.table(ddsaSFNormcounts_chr_gene1, "ddsaSFNormcounts_chr_gene1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)

#All  genes mat ratio
ddsaSFNormcounts_chr_gene1_ratio <- ddsaSFNormcounts_chr_gene1
ddsaSFNormcounts_chr_gene1_ratio["WT_allele_ratio"] <- ddsaSFNormcounts_chr_gene1_ratio$JB1_WT_JF1 / (ddsaSFNormcounts_chr_gene1_ratio$JB1_WT_B6 + ddsaSFNormcounts_chr_gene1_ratio$JB1_WT_JF1)
ddsaSFNormcounts_chr_gene1_ratio["ZFP57KO_allele_ratio"] <- ddsaSFNormcounts_chr_gene1_ratio$JB1_ZFP57KO_JF1 / (ddsaSFNormcounts_chr_gene1_ratio$JB1_ZFP57KO_B6 + ddsaSFNormcounts_chr_gene1_ratio$JB1_ZFP57KO_JF1)
head(ddsaSFNormcounts_chr_gene1_ratio)
write.table(ddsaSFNormcounts_chr_gene1_ratio, "ddsaSFNormcounts_chr_gene1_ratio.txt", sep = "\t", quote = F, append = F, row.names = F)



#I created a same sheet: 
fgrep -f basilia_genes_for_dotplot.txt ddsaSFNormcounts_chr_gene1.txt -w | sort -k4,4 > basilia_genes_for_dotplot_ddsaSFNormcounts.txt
#Only 29 genes covered in T0 data and 6 genes (A230057D06Rik, Peg12, Airn, Igf2os, Dcn, Th) has no expression
#Manually done
#Open basilia_genes_for_dotplot_ddsaSFNormcounts.txt and if 0 occurs in both allele, we replaced it with 1 (atleast one read), the other allele was still high
cp basilia_genes_for_dotplot_ddsaSFNormcounts.txt basilia_genes_for_dotplot_ddsaSFNormcounts_re.txt
basilia_genes_for_dotplot_ddsaSFNormcounts <- read.table("basilia_genes_for_dotplot_ddsaSFNormcounts_re.txt", header = F)
colnames(basilia_genes_for_dotplot_ddsaSFNormcounts) <- colnames(ddsaSFNormcounts_chr_gene1)
head(basilia_genes_for_dotplot_ddsaSFNormcounts)
basilia_genes_for_dotplot_ddsaSFNormcounts["WT_allele_ratio"] <- basilia_genes_for_dotplot_ddsaSFNormcounts$JB1_WT_JF1 / (basilia_genes_for_dotplot_ddsaSFNormcounts$JB1_WT_B6 + basilia_genes_for_dotplot_ddsaSFNormcounts$JB1_WT_JF1)
basilia_genes_for_dotplot_ddsaSFNormcounts["ZFP57KO_allele_ratio"] <- basilia_genes_for_dotplot_ddsaSFNormcounts$JB1_ZFP57KO_JF1 / (basilia_genes_for_dotplot_ddsaSFNormcounts$JB1_ZFP57KO_B6 + basilia_genes_for_dotplot_ddsaSFNormcounts$JB1_ZFP57KO_JF1)

data_imp_gene_indiv_asp_norm <- basilia_genes_for_dotplot_ddsaSFNormcounts
# order
data_imp_gene_indiv_asp_norm <- data_imp_gene_indiv_asp_norm[order(data_imp_gene_indiv_asp_norm$WT_allele_ratio),]
head(data_imp_gene_indiv_asp_norm)
dim(data_imp_gene_indiv_asp_norm)
write.table(data_imp_gene_indiv_asp_norm, "data_imp_gene_indiv_asp_norm.txt", sep = "\t", quote = F, append = F, row.names = F)
head(data_imp_gene_indiv_asp_norm,1)
summary(data_imp_gene_indiv_asp_norm[,c(10,11)])
data_imp_gene_indiv_asp_norm_re <- data_imp_gene_indiv_asp_norm[,c(4,10,11)]
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
ggsave("barplot_data_imp_gene_indiv_asp_norm_rebarplot1_half.svg", width=13, height=20, units="cm", dpi=96)
ggsave("barplot_data_imp_gene_indiv_asp_norm_rebarplot1_half.jpg", width=13, height=20, units="cm", dpi=96)



pdot <-ggplot(data=data_imp_gene_indiv_asp_norm_rebarplot1_half, aes(x=Gene, y=Ratio, label=Group, color="black", fill=Group)) + 
  geom_hline(yintercept = c(-0.17,0,0.19), colour = "grey", linetype=c("dashed","solid","dashed")) + 
  geom_dotplot(binaxis='y', stackdir='centerwhole', stackratio=1.5, dotsize=1.2) + theme_classic() 
pdot + coord_flip()+ ylim(c(-0.5,0.5))+
  scale_fill_manual(values=c("green","red"))+
  scale_color_manual(values=c("black","black","black","black")) 

ggsave("dotplot_data_imp_gene_indiv_asp_norm_sortrebarplot1_half.svg", width=13, height=20, units="cm", dpi=96)
ggsave("dotplot_data_imp_gene_indiv_asp_norm_sortrebarplot1_half.jpg", width=13, height=20, units="cm", dpi=96)

pdotsquare <-ggplot(data=data_imp_gene_indiv_asp_norm_rebarplot1_half, aes(x=Gene, y=Ratio, label=Group, color="black", fill=Group)) + 
  geom_hline(yintercept = c(-0.17,0,0.19), colour = "grey", linetype=c("dashed","solid","dashed")) + 
  geom_point(aes(shape=Group, color=Group, fill=Group, size=Group), position = position_dodge(width = 0.50)) + theme_classic()  

pdotsquare + coord_flip()+ ylim(c(-0.5,0.5))+
  scale_fill_manual(values=c("green","red"))+
  scale_color_manual(values=c("black","black","black","black")) +
  scale_shape_manual(values=c(22,22))+
  scale_size_manual(values=c(1.5,1.5,1.5,1.5)) 

ggsave("trisquareplot_data_imp_gene_indiv_asp_norm_sortrebarplot1_halfdodge.svg", width=13, height=20, units="cm", dpi=96)
ggsave("trisquareplot_data_imp_gene_indiv_asp_norm_sortrebarplot1_halfdodge.jpg", width=13, height=20, units="cm", dpi=96)

pdotsquare <-ggplot(data=data_imp_gene_indiv_asp_norm_rebarplot1_half, aes(x=Gene, y=Ratio, label=Group, color="black", fill=Group)) + 
  geom_hline(yintercept = c(-0.17,0,0.19), colour = "grey", linetype=c("dashed","solid","dashed")) + 
  geom_point(aes(shape=Group, color=Group, fill=Group, size=Group), position = position_dodge(width = 0.1)) + theme_classic()  

pdotsquare + coord_flip()+ ylim(c(-0.5,0.5))+
  scale_fill_manual(values=c("green","red"))+
  scale_color_manual(values=c("black","black","black","black")) +
  scale_shape_manual(values=c(22,22))+
  scale_size_manual(values=c(1.5,1.5,1.5,1.5)) 

ggsave("trisquareplot_data_imp_gene_indiv_asp_norm_sortrebarplot1_half.svg", width=13, height=20, units="cm", dpi=96)
ggsave("trisquareplot_data_imp_gene_indiv_asp_norm_sortrebarplot1_half.jpg", width=13, height=20, units="cm", dpi=96)


#Fisher test and tag chromsome positions
#ddsaSFNormcounts_chr_gene1 <- read.table("ddsaSFNormcounts_chr_gene1.txt", header = T)
head(ddsaSFNormcounts_chr_gene1)
dim(ddsaSFNormcounts_chr_gene1)
rownames(ddsaSFNormcounts_chr_gene1) <- paste(ddsaSFNormcounts_chr_gene1$Genes,ddsaSFNormcounts_chr_gene1$ensid, sep = "%")
head(ddsaSFNormcounts_chr_gene1)
ddsaSFNormcounts_chr_gene1_fisher <- ddsaSFNormcounts_chr_gene1[,6:9]
head(ddsaSFNormcounts_chr_gene1_fisher)
ddsaSFNormcounts_chr_gene1_fisher1 <- apply(ddsaSFNormcounts_chr_gene1_fisher,1, function(x) fisher.test(matrix(x,nrow =2))$p.value)
#warning will be observed because we need only integers for runnung this function.
#In fisher.test(matrix(x, nrow = 2)) :'x' has been rounded to integer: Mean relative difference:
head(ddsaSFNormcounts_chr_gene1_fisher1)
ddsaSFNormcounts_chr_gene1_fisher["fisherpvalue"] <- data.frame(ddsaSFNormcounts_chr_gene1_fisher1)
head(ddsaSFNormcounts_chr_gene1_fisher)
ddsaSFNormcounts_chr_gene1_fisher["fisheradjpval"] <- p.adjust(ddsaSFNormcounts_chr_gene1_fisher1,method="BH")
head(ddsaSFNormcounts_chr_gene1_fisher)
dim(ddsaSFNormcounts_chr_gene1_fisher)
write.table(ddsaSFNormcounts_chr_gene1_fisher, "ddsaSFNormcounts_chr_gene1_fisher.txt", sep="\t", quote = FALSE, append = FALSE)

ddsaSFNormcounts_chr_gene1_fishercoordinate <- cbind.data.frame(ddsaSFNormcounts_chr_gene1, ddsaSFNormcounts_chr_gene1_fisher)
head(ddsaSFNormcounts_chr_gene1_fishercoordinate)
write.table(ddsaSFNormcounts_chr_gene1_fishercoordinate, "ddsaSFNormcounts_chr_gene1_fishercoordinate.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)
writexl::write_xlsx(ddsaSFNormcounts_chr_gene1_fishercoordinate, "JB1_ESCs_T0_ddsaSFNormcounts_chr_gene1_fishercoordinate_mm10.xlsx")


#WithID
ddsaSFNormcounts_chr_gene1_fisher["GeneID"] <- rownames(ddsaSFNormcounts_chr_gene1_fisher)
ddsaSFNormcounts_chr_gene1_fisher_0.05 <- ddsaSFNormcounts_chr_gene1_fisher[which(ddsaSFNormcounts_chr_gene1_fisher$fisheradjpval <0.05),]
head(ddsaSFNormcounts_chr_gene1_fisher_0.05)
dim(ddsaSFNormcounts_chr_gene1_fisher_0.05)
write.table(ddsaSFNormcounts_chr_gene1_fisher_0.05, "ddsaSFNormcounts_chr_gene1_fisher_0.05.txt", sep="\t", quote = FALSE, append = FALSE)

#WithID
ddsaSFNormcounts_chr_gene1_fisher_0.1 <- ddsaSFNormcounts_chr_gene1_fisher[which(ddsaSFNormcounts_chr_gene1_fisher$fisheradjpval <0.1),]
head(ddsaSFNormcounts_chr_gene1_fisher_0.1)
dim(ddsaSFNormcounts_chr_gene1_fisher_0.1)
write.table(ddsaSFNormcounts_chr_gene1_fisher_0.1, "ddsaSFNormcounts_chr_gene1_fisher_0.1.txt", sep="\t", quote = FALSE, append = FALSE)

#After Prof's suggestion only the statistically significant genes (fisher/BH adj p-value  < 0.1) were selected 
basilia_genes_for_dotplot <- read.table("basilia_genes_for_dotplot.txt")
colnames(basilia_genes_for_dotplot) <- c("Genes")
basilia_imp_genes_ddsSFNormcounts_chr.pos = merge(basilia_genes_for_dotplot, ddsaSFNormcounts_chr_gene1_fishercoordinate, by="Genes", all.x=FALSE)
head(basilia_imp_genes_ddsSFNormcounts_chr.pos)
dim(basilia_imp_genes_ddsSFNormcounts_chr.pos)

#Select Fisher/BH adj p-value < 0.05
basilia_imp_genes_ddsSFNormcounts_chr_sig <- basilia_imp_genes_ddsSFNormcounts_chr.pos[which(basilia_imp_genes_ddsSFNormcounts_chr.pos$fisheradjpval < 0.05),]
dim(basilia_imp_genes_ddsSFNormcounts_chr_sig)
head(basilia_imp_genes_ddsSFNormcounts_chr_sig)
#Also remove Cdkn1c and Kcnq1ot1
basilia_imp_genes_ddsSFNormcounts_chr_resig <- basilia_imp_genes_ddsSFNormcounts_chr_sig[which(basilia_imp_genes_ddsSFNormcounts_chr_sig$Genes != "Cdkn1c"),] 
basilia_imp_genes_ddsSFNormcounts_chr_resig <- basilia_imp_genes_ddsSFNormcounts_chr_resig[which(basilia_imp_genes_ddsSFNormcounts_chr_resig$Genes != "Kcnq1ot1"),]
dim(basilia_imp_genes_ddsSFNormcounts_chr_resig) #13 left


basilia_imp_genes_ddsSFNormcounts_chr_resig["WT_allele_ratio"] <- basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_WT_JF1 / (basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_WT_B6 + basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_WT_JF1)
basilia_imp_genes_ddsSFNormcounts_chr_resig["ZFP57KO_allele_ratio"] <- basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_ZFP57KO_JF1 / (basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_ZFP57KO_B6 + basilia_imp_genes_ddsSFNormcounts_chr_resig$JB1_ZFP57KO_JF1)

data_imp_gene_indiv_asp_norm_resig <- basilia_imp_genes_ddsSFNormcounts_chr_resig
# order
data_imp_gene_indiv_asp_norm_resig <- data_imp_gene_indiv_asp_norm_resig[order(data_imp_gene_indiv_asp_norm_resig$WT_allele_ratio),]
head(data_imp_gene_indiv_asp_norm_resig)
dim(data_imp_gene_indiv_asp_norm_resig)
write.table(data_imp_gene_indiv_asp_norm_resig, "data_imp_gene_indiv_asp_norm_resig.txt", sep = "\t", quote = F, append = F, row.names = F)
head(data_imp_gene_indiv_asp_norm_resig,1)
summary(data_imp_gene_indiv_asp_norm_resig[,c(16,17)])
data_imp_gene_indiv_asp_norm_resig_re <- data_imp_gene_indiv_asp_norm_resig[,c(1,16,17)]
head(data_imp_gene_indiv_asp_norm_resig_re)
colnames(data_imp_gene_indiv_asp_norm_resig_re) <- c("Gene", "WT_avg_allele_ratio", "ZFP57KO_avg_allele_ratio")
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
  geom_hline(yintercept = c(-0.17,0,0.19), colour = "grey", linetype=c("dashed","solid","dashed")) + 
  geom_point(aes(shape=Group, color=Group, fill=Group, size=Group), position = position_dodge(width = 0.1)) + theme_classic()  

psigdotsquare + coord_flip()+ ylim(c(-0.5,0.5))+
  scale_fill_manual(values=c("#3CB371","#DC143C"))+
  scale_color_manual(values=c("black","black")) +
  scale_shape_manual(values=c(21,21))+
  scale_size_manual(values=c(2,2)) 

ggsave("trisquareplot_data_imp_gene_indiv_asp_norm_resig_sortrebarplot1_half.svg", width=13, height=15, units="cm", dpi=96)
ggsave("trisquareplot_data_imp_gene_indiv_asp_norm_resig_sortrebarplot1_half.jpg", width=13, height=15, units="cm", dpi=96)



#Prop test and tag chromsome positions
ddsaSFNormcounts_chr_gene1 <- ddsSFNormcounts_chr.pos
head(ddsaSFNormcounts_chr_gene1)
library(metaseqR)
head(ddsaSFNormcounts_chr_gene1)
dim(ddsaSFNormcounts_chr_gene1)
ddsaSFNormcounts_chr_gene1["id"] <- rownames(ddsaSFNormcounts_chr_gene1)
rownames(ddsaSFNormcounts_chr_gene1) <- paste(ddsaSFNormcounts_chr_gene1$Genes,ddsaSFNormcounts_chr_gene1$ensid, sep = "%")
ddsaSFNormcounts_chr_gene1_prop <- ddsaSFNormcounts_chr_gene1[,6:9]
head(ddsaSFNormcounts_chr_gene1_prop)
ddsaSFNormcounts_chr_gene1_prop["JB1_WT_alavg"] <- ddsaSFNormcounts_chr_gene1_prop$JB1_WT_B6 + ddsaSFNormcounts_chr_gene1_prop$JB1_WT_JF1
ddsaSFNormcounts_chr_gene1_prop["JB1_ZFP57KO_alavg"] <- ddsaSFNormcounts_chr_gene1_prop$JB1_ZFP57KO_B6 + ddsaSFNormcounts_chr_gene1_prop$JB1_ZFP57KO_JF1

#Take allelic ratio
ddsaSFNormcounts_chr_gene1_prop["WT_allele_ratio"] <- ddsaSFNormcounts_chr_gene1_prop$JB1_WT_JF1 / (ddsaSFNormcounts_chr_gene1_prop$JB1_WT_B6 + ddsaSFNormcounts_chr_gene1_prop$JB1_WT_JF1)
ddsaSFNormcounts_chr_gene1_prop["ZFP57KO_allele_ratio"] <- ddsaSFNormcounts_chr_gene1_prop$JB1_ZFP57KO_JF1 / (ddsaSFNormcounts_chr_gene1_prop$JB1_ZFP57KO_B6 + ddsaSFNormcounts_chr_gene1_prop$JB1_ZFP57KO_JF1)

head(ddsaSFNormcounts_chr_gene1_prop)
#Filter B6+JF1 > 10 in WT, if both the alleles are 10 in either of the replicates
ddsaSFNormcounts_chr_gene1_prop.filt <- ddsaSFNormcounts_chr_gene1_prop[which(ddsaSFNormcounts_chr_gene1_prop$JB1_WT_alavg > 10 &
                                                                                ddsaSFNormcounts_chr_gene1_prop$JB1_ZFP57KO_alavg > 0),]


dim(ddsaSFNormcounts_chr_gene1_prop.filt)
head(ddsaSFNormcounts_chr_gene1_prop.filt,1)
summary(ddsaSFNormcounts_chr_gene1_prop.filt)
#WT
JB1_WT_alavg.prop <- Map(prop.test,x =ddsaSFNormcounts_chr_gene1_prop.filt$JB1_WT_B6, n= ddsaSFNormcounts_chr_gene1_prop.filt$JB1_WT_alavg, p=0.5)

ddsaSFNormcounts_chr_gene1_prop.filt["JB1_WT_alavg.pvalue"] <- data.frame(capture.output(for (i in seq_along(JB1_WT_alavg.prop)){
  cat(JB1_WT_alavg.prop[[i]]$p.value, "\n")
}))


#No need to Combine p-values calculated from prop.test 

head(ddsaSFNormcounts_chr_gene1_prop.filt)
dim(ddsaSFNormcounts_chr_gene1_prop.filt)
tail(ddsaSFNormcounts_chr_gene1_prop.filt)
ddsaSFNormcounts_chr_gene1_prop.filt <- data.frame(ddsaSFNormcounts_chr_gene1_prop.filt)
#Apply BH on p-value from prop.test
#Important: The p-value obtained are factors and need to be converted to numeric as follows as.numeric(as.character(p-value)
ddsaSFNormcounts_chr_gene1_prop.filt["JB1_WT_alavg.fdr"] <- p.adjust(as.numeric(as.character(ddsaSFNormcounts_chr_gene1_prop.filt$JB1_WT_alavg.pvalue)), method="BH")
head(ddsaSFNormcounts_chr_gene1_prop.filt)
dim(ddsaSFNormcounts_chr_gene1_prop.filt)
ddsaSFNormcounts_chr_gene1_prop.filt_wt <- ddsaSFNormcounts_chr_gene1_prop.filt
ddsaSFNormcounts_chr_gene1_prop.filt_wt["id"] <- rownames(ddsaSFNormcounts_chr_gene1_prop.filt_wt)
writexl::write_xlsx(ddsaSFNormcounts_chr_gene1_prop.filt_wt, "ddsaSFNormcounts_chr_gene1_prop.filt_wt.xlsx")
write.table(ddsaSFNormcounts_chr_gene1_prop.filt_wt, "ddsaSFNormcounts_chr_gene1_prop.filt_wt.txt", sep = "\t", quote = F, append = F)
#manually also checked

#ZFP57_KO
JB1_ZFP57KO_alavg.prop <- Map(prop.test,x =ddsaSFNormcounts_chr_gene1_prop.filt$JB1_ZFP57KO_B6, n= ddsaSFNormcounts_chr_gene1_prop.filt$JB1_ZFP57KO_alavg, p=0.5)

ddsaSFNormcounts_chr_gene1_prop.filt["JB1_ZFP57KO_alavg.pvalue"] <- data.frame(capture.output(for (i in seq_along(JB1_ZFP57KO_alavg.prop)){
  cat(JB1_ZFP57KO_alavg.prop[[i]]$p.value, "\n")
}))

#No need to Combine p-values calculated from prop.test 

head(ddsaSFNormcounts_chr_gene1_prop.filt)
tail(ddsaSFNormcounts_chr_gene1_prop.filt)
#Apply BH on p-value from prop.test
ddsaSFNormcounts_chr_gene1_prop.filt["JB1_ZFP57KO_alavg.fdr"] <- p.adjust(as.numeric(as.character(ddsaSFNormcounts_chr_gene1_prop.filt$JB1_ZFP57KO_alavg.pvalue)),method="BH")
head(ddsaSFNormcounts_chr_gene1_prop.filt)


#Introduce gene name
ddsaSFNormcounts_chr_gene1_prop.filt["id"] <- data.frame(rownames(ddsaSFNormcounts_chr_gene1_prop.filt))
head(ddsaSFNormcounts_chr_gene1_prop.filt)
library(splitstackshape)
ddsaSFNormcounts_chr_gene1_prop.filt.sep <- cSplit(ddsaSFNormcounts_chr_gene1_prop.filt, "id", "%")
dim(ddsaSFNormcounts_chr_gene1_prop.filt.sep)
head(ddsaSFNormcounts_chr_gene1_prop.filt.sep)
ddsaSFNormcounts_chr_gene1_prop.filt_gene = data.frame(ddsaSFNormcounts_chr_gene1_prop.filt.sep[,c(13,14,1:12)])
head(ddsaSFNormcounts_chr_gene1_prop.filt_gene)
write.table(ddsaSFNormcounts_chr_gene1_prop.filt_gene, "ddsaSFNormcounts_chr_gene1_prop.filt_gene.txt", sep = "\t", quote = F, row.names = F, append = F)
writexl::write_xlsx(ddsaSFNormcounts_chr_gene1_prop.filt_gene,"AlleleSpecific_ddsaSFNormcounts_chr_gene1_proptest_all.xlsx")

#Assign chr to ddsaSFNormcounts_chr_gene1_prop.filt_gene
chr.position  =  read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt",header=FALSE)
head(chr.position)
colnames(chr.position) <- c("id_2", "Genes", "chr", "start", "end")
head(chr.position)
dim(chr.position)
head(ddsaSFNormcounts_chr_gene1_prop.filt_gene)
dim(ddsaSFNormcounts_chr_gene1_prop.filt_gene)
ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr = merge(ddsaSFNormcounts_chr_gene1_prop.filt_gene, chr.position, by="id_2", all.x=FALSE)
head(ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr)
dim(ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr)
ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr <- ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr[,c(16:18,1:15)]
write.table(ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr, "ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr.txt", sep = "\t", quote = F, row.names = F, append = F)



#Monoalleically expressed in Wildtype
ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig <- ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr[which(ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr$JB1_WT_alavg.fdr<0.05),]
head(ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig)
dim(ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig)
write.table(ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig, "ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig.txt", sep = "\t", quote = F, row.names = F, append = F)
writexl::write_xlsx(ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig,"AlleleSpecific_ddsaSFNormcounts_chr_gene1_proptest_WT_fdrsig.xlsx")
ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig0.1 <- ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr[which(ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr$JB1_WT_alavg.fdr<0.1),]
head(ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig0.1)
dim(ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig0.1)
write.table(ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig0.1, "ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig0.1.txt", sep = "\t", quote = F, row.names = F, append = F)
writexl::write_xlsx(ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig0.1,"AlleleSpecific_ddsaSFNormcounts_chr_gene1_proptest_WT_fdrsig0.1.xlsx")


#Monoalleically expressed in ZFP5KO
ddsaSFNormcounts_chr_gene1_prop.filt_gene_ZFP57KO_fdrsig <- ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr[which(ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr$JB1_ZFP57KO_alavg.fdr<0.05),]
head(ddsaSFNormcounts_chr_gene1_prop.filt_gene_ZFP57KO_fdrsig)
dim(ddsaSFNormcounts_chr_gene1_prop.filt_gene_ZFP57KO_fdrsig)
write.table(ddsaSFNormcounts_chr_gene1_prop.filt_gene_ZFP57KO_fdrsig, "ddsaSFNormcounts_chr_gene1_prop.filt_gene_ZFP57KO_fdrsig.txt", sep = "\t", quote = F, row.names = F, append = F)
writexl::write_xlsx(ddsaSFNormcounts_chr_gene1_prop.filt_gene_ZFP57KO_fdrsig,"AlleleSpecific_ddsaSFNormcounts_chr_gene1_proptest_ZFP57KO_fdrsig.xlsx")
ddsaSFNormcounts_chr_gene1_prop.filt_gene_ZFP57KO_fdrsig0.1 <- ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr[which(ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr$JB1_ZFP57KO_alavg.fdr<0.1),]
head(ddsaSFNormcounts_chr_gene1_prop.filt_gene_ZFP57KO_fdrsig0.1)
dim(ddsaSFNormcounts_chr_gene1_prop.filt_gene_ZFP57KO_fdrsig0.1)
write.table(ddsaSFNormcounts_chr_gene1_prop.filt_gene_ZFP57KO_fdrsig0.1, "ddsaSFNormcounts_chr_gene1_prop.filt_gene_ZFP57KO_fdrsig0.1.txt", sep = "\t", quote = F, row.names = F, append = F)
writexl::write_xlsx(ddsaSFNormcounts_chr_gene1_prop.filt_gene_ZFP57KO_fdrsig0.1,"AlleleSpecific_ddsaSFNormcounts_chr_gene1_proptest_ZFP57KO_fdrsig0.1.xlsx")

#Assign Imprinted gene names
fgrep -f imprinted_gene_name.txt ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig.txt -w > imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig.txt 
#remove chrX, (chrY, chrM) already not present in  imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig.txt
grep chrX imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig.txt -v > imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX.txt
imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX <- read.table("imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX.txt", header = F, stringsAsFactors = F)
head(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX)
colnames(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX) <- colnames(ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig)
head(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX)
dim(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX)
writexl::write_xlsx(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX,"imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX.xlsx")

#Add id as it is required for fisher merge
imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX["id"] <- imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX$id_2
dim(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX)
head(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX)
#Add ZFP57/KAP1 peaks  information
head(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re)
imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher <- merge(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re, imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX, by="id", all.x=F)
head(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher,1)
dim(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher)
writexl::write_xlsx(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher,"imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher.xlsx")

#Note in this data.frame fisherpvalue fisheradjpval are allelic disbalance 2x2 matrix fisher test while *.comb.fisherpvalue are from prop test p-value combine, They are different
#Now all imprinted genes below are statisitcal significant from Prop test in WT
#Filter by allelic ratio
imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af <- imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher[which(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher$WT_allele_ratio <= 0.33 |
                                                                                                                                                                                                                                    imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher$WT_allele_ratio >= 0.67),]
head(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af)
dim(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af)
writexl::write_xlsx(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af,"ESC_imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher.xlsx")

imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af2 <- imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af[which(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af$Genes!="Nnat"),]
head(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af2)
dim(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af2)

#Re-Dotplot, only significant prop test and sufficiently expressed,  imprinted genes
data_for_merge_dotplot_ESCs <- imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af2[,c(6,34,35,37)]
head(data_for_merge_dotplot_ESCs)
dim(data_for_merge_dotplot_ESCs)
colnames(data_for_merge_dotplot_ESCs) <- c("ESCs_Genes", "ESCs_WT_allele_ratio", "ESCs_ZFP57KO_allele_ratio", "ESCs_JB1_WT_alavg.comb.fdr")
data_for_merge_dotplot_ESCs["Genes"] <- data_for_merge_dotplot_ESCs$ESCs_Genes
write.table(data_for_merge_dotplot_ESCs, "data_for_merge_dotplot_ESCs.txt", sep = "\t", quote = F, row.names = F, append = F)
nonsig_genesdata_for_merge_dotplot_ESCs <- data_for_merge_dotplot_ESCs[which(data_for_merge_dotplot_ESCs$ESCs_JB1_WT_alavg.comb.fdr > 0.05),]
write.table(nonsig_genesdata_for_merge_dotplot_ESCs$Genes, "nonsig_genesdata_for_merge_dotplot_ESCs.geneid", sep = "\t", quote = F, row.names = F, append = F)
dim(data_for_merge_dotplot_ESCs)


#Assign chr to ddsaSFNormcounts_chr_gene1_prop.filt_gene
chr.position  =  read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt",header=FALSE)
head(chr.position)
colnames(chr.position) <- c("id_3", "Genes", "chr", "start", "end")
head(chr.position)
dim(chr.position)
head(ddsaSFNormcounts_chr_gene1_prop.filt_gene)
dim(ddsaSFNormcounts_chr_gene1_prop.filt_gene)
ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr = merge(ddsaSFNormcounts_chr_gene1_prop.filt_gene, chr.position, by="id_3", all.x=FALSE)
head(ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr)
dim(ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr)
ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr <- ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr[,c(16:18,1:14)]


#Do it for all prop genes
#Add id as it is required for fisher merge
ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr["id"] <- ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr$id_2
dim(ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr)
head(ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr)

#Add strand information
#Add closest distance to Deseq2 noramlized counts
grep ENS ddsaSFNormcounts_chr_gene1.txt | sort -k1,1 -k2,2n > ddsaSFNormcounts_chr_gene1.sort.txt
#While intersecting with bedtools remove duplication which happen due to overlap of Zfp57 peaks, keep only once. It is even Ok for imprinted genes even if they overlap more than one time with peaks
bedtools closest -a ddsaSFNormcounts_chr_gene1.sort.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | sort -k5,5 -u > ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10.txt
#Import this file and perform Fisher exact test and then combine
ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10 <- read.table("ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10.txt", header = F)
head(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10)
dim(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10) 
colnames(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10) <- c("chr", "start", "end", "Genes","id","JB1_WT_B6","JB1_WT_JF1","JB1_ZFP57_KO_B6","JB1_ZFP57_KO_JF1","chrPeak","startpeak","endPeak","peakID","qualityPeak","ClosestDistance")
head(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10)
rownames(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10) <- paste(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10$Genes,ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10$id, sep = "%")
partialgddsSFNormcounts <- ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10
#Take average
partialddsSFNormcountsavg <- partialgddsSFNormcounts[,c(4,6:9)]
head(partialddsSFNormcountsavg)
colnames(partialddsSFNormcountsavg) <- c("GeneID","WT.B6","WT.JF1","ZFP57KO.B6","ZFP57KO.JF1")
#Fisher exact Test
partialddsSFNormcountsavg_fisher <- partialddsSFNormcountsavg[,2:5]
head(partialddsSFNormcountsavg_fisher)
partialddsSFNormcountsavg_fisher["fisherpvalue"] <- apply(partialddsSFNormcountsavg_fisher,1, function(x) fisher.test(matrix(x,nrow =2))$p.value)

rownames(partialddsSFNormcountsavg_fisher) <- rownames(partialddsSFNormcountsavg)
partialddsSFNormcountsavg_fisher["GeneID"] <- partialddsSFNormcountsavg$GeneID
head(partialddsSFNormcountsavg_fisher)
#warning will be observed because we need only integers for runnung this function.
#In fisher.test(matrix(x, nrow = 2)) :'x' has been rounded to integer: Mean relative difference:

#Recombine full matrix
ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher <- cbind.data.frame(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10,partialddsSFNormcountsavg_fisher)
head(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher,1)
ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher["fisheradjpval"] <- p.adjust(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher$fisherpvalue,method="BH")
head(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher)
write.table(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher, "ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.txt", sep="\t", quote = FALSE, append = FALSE)
#Filter p<0.05
ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher_0.05 <- ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher[which(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher$fisheradjpval <0.05),]
head(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher_0.05)
dim(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher_0.05)
write.table(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher_0.05, "ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher_0.05.txt", sep="\t", quote = FALSE, append = FALSE)

ens_gene_names_chrpos_dedup_M20.strand <- read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.strand.txt")
head(ens_gene_names_chrpos_dedup_M20.strand)
colnames(ens_gene_names_chrpos_dedup_M20.strand) <- c("chr", "start", "end","gene","strand","id","Gene")
dim(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher)
head(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher)
ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes = merge(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher, ens_gene_names_chrpos_dedup_M20.strand, by="id", all.x=TRUE)
head(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes)
dim(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes)
ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re <- ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes[,c(23:25,27:28,1,6:20,22)]
head(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re)
dim(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re)
writexl::write_xlsx(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re, "AlleleSpecific_gddsaSFNormcounts_overlapped_ZFP57_KAP1_peaks_fisher.genes.re.xlsx")

ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr$id <- ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr$id_2

head(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re)
ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher <- merge(ddsaSFNormcounts_chr_gene1.sort_JB1_Zfp57_overlapped_Kap1_mm10_fisher.genes.re, ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr, by="id", all.x=TRUE)
head(ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher,1)
dim(ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher)
writexl::write_xlsx(ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher,"ddsaSFNormcounts_chr_gene1_prop.filt_gene_chr_JB1_Zfp57_overlapped_Kap1_mm10_fisher.xlsx")
#Gene with NA in prop data column are not covered in fisher test, as checked manually


#Search npc sig genes in ESC created a same sheet: 
fgrep -f data_for_merge_dotplot_NPCs_Genes_sigsel.txt ddsaSFNormcounts_chr_gene1.txt -w  > selectednpcsig_genes_for_dotplot_ddsaSFNormcounts.txt
#Only 29 genes covered in T0 data and 3 genes (A230057D06Rik, Peg12, Igf2os) has no expression in ESCs (ND in table)

selectednpcsig_genes_for_dotplot_ddsaSFNormcounts <- read.table("selectednpcsig_genes_for_dotplot_ddsaSFNormcounts.txt", header = F)
head(selectednpcsig_genes_for_dotplot_ddsaSFNormcounts)
colnames(selectednpcsig_genes_for_dotplot_ddsaSFNormcounts) <- colnames(ddsaSFNormcounts_chr_gene1[,1:(length(ddsaSFNormcounts_chr_gene1)-1)])
head(selectednpcsig_genes_for_dotplot_ddsaSFNormcounts)
selectednpcsig_genes_for_dotplot_ddsaSFNormcounts["WT_allele_ratio"] <- selectednpcsig_genes_for_dotplot_ddsaSFNormcounts$JB1_WT_JF1 / (selectednpcsig_genes_for_dotplot_ddsaSFNormcounts$JB1_WT_B6 + selectednpcsig_genes_for_dotplot_ddsaSFNormcounts$JB1_WT_JF1)
selectednpcsig_genes_for_dotplot_ddsaSFNormcounts["ZFP57KO_allele_ratio"] <- selectednpcsig_genes_for_dotplot_ddsaSFNormcounts$JB1_ZFP57KO_JF1 / (selectednpcsig_genes_for_dotplot_ddsaSFNormcounts$JB1_ZFP57KO_B6 + selectednpcsig_genes_for_dotplot_ddsaSFNormcounts$JB1_ZFP57KO_JF1)
head(selectednpcsig_genes_for_dotplot_ddsaSFNormcounts)
dim(selectednpcsig_genes_for_dotplot_ddsaSFNormcounts)
write.table(selectednpcsig_genes_for_dotplot_ddsaSFNormcounts, "selectednpcsig_genes_for_dotplot_ddsaSFNormcounts1.txt", sep = "\t", quote = F, row.names = F, append = F)
#remove 10 low expressed genes Nnat, Sgce, Asb4, Ndn, Magel2, Mkrn3, Zim1, Igf2, Cdkn1c, Dlk1 (LE in table), See UCSC visualization selection https://docs.google.com/presentation/d/1urcxQbkGZNhl9k8wVDpmKPY-P2kDpAVjm00P_QCCk4g/edit#slide=id.gbb3bb6bd07_3_32
fgrep -f rmv_genes_for_ESCs.txt selectednpcsig_genes_for_dotplot_ddsaSFNormcounts1.txt -v > selectednpcsig_genes_for_dotplot_ddsaSFNormcounts1_rmv.txt
#19 imprinted genes left in ESC
selectednpcsig_dotplot_ESCs <- read.table("selectednpcsig_genes_for_dotplot_ddsaSFNormcounts1_rmv.txt", header = T, stringsAsFactors = F)
head(selectednpcsig_dotplot_ESCs)
dim(selectednpcsig_dotplot_ESCs)
selectednpcsig_dotplot_ESCs <- selectednpcsig_dotplot_ESCs[,c(4,10,11)]
head(selectednpcsig_dotplot_ESCs)
colnames(selectednpcsig_dotplot_ESCs) <- c("ESCs_Genes", "ESCs_WT_allele_ratio", "ESCs_ZFP57KO_allele_ratio")
selectednpcsig_dotplot_ESCs["Genes"] <- selectednpcsig_dotplot_ESCs$ESCs_Genes
head(selectednpcsig_dotplot_ESCs)
dim(selectednpcsig_dotplot_ESCs)
write.table(selectednpcsig_dotplot_ESCs, "selectednpcsig_dotplot_ESCs.txt", sep = "\t", quote = F, row.names = F, append = F)

#So first observation is statistically significant genes from prop test in ESCs
#Significant_imprinted_genes and  Remove chrX genes (male cell line)
dim(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX)#42

#Allelic Ratio filter 0.33=<AR>=0.67
dim(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af) #30

#Remove very  low expressed genes ==Nnat
dim(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af2)#29
#or 
dim(data_for_merge_dotplot_ESCs)#29

#Overlap this 29 genes with NPC sig genes and select genes which were significant in NPC
data_for_merge_dotplot_NPCs_Genes_sigsel <- read.table("data_for_merge_dotplot_NPCs_Genes_sigsel.txt", header = F)
colnames(data_for_merge_dotplot_NPCs_Genes_sigsel) <- "Genes"
imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af3 <- merge(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af2,data_for_merge_dotplot_NPCs_Genes_sigsel, by="Genes", all.x=FALSE, all.y=FALSE)
dim(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af3)
writexl::write_xlsx(imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af3,
                    "imprinted_ddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af3.xlsx")

#Out of 29  statistically signifcant imrinted genes, and after ratio filtering and low count filtering in ESC (prop test) only 19
#were matched with NPC selected genes. 
#The other 10 were not present in NPC selected genes: Tspan32, Cobl, Igf2r, Ano1, Magi2, Dhcr7, Art5, Zim2, Zim3, Zfp264.



############################################# END OF ANALYSIS ###############################################3
#SNPs inside exons
gencode.vM20.chr_patch_hapl_scaff.annotation.gtf <- rtracklayer::import('/home/ankitv/ref_av/gencodes/gencode_M20/prep/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf')
gencode.vM20.chr_patch_hapl_scaff.annotation.gtf_df =as.data.frame(gencode.vM20.chr_patch_hapl_scaff.annotation.gtf)
head(gencode.vM20.chr_patch_hapl_scaff.annotation.gtf_df)
dim(gencode.vM20.chr_patch_hapl_scaff.annotation.gtf_df)
exons.vM20.annotation <- gencode.vM20.chr_patch_hapl_scaff.annotation.gtf_df[which(gencode.vM20.chr_patch_hapl_scaff.annotation.gtf_df$type == "exon"),]
head(exons.vM20.annotation)
dim(exons.vM20.annotation)

#Count exons
exons.vM20.annotation_gene_chr <- exons.vM20.annotation[,c(1,2,3,5,7,10,12,15,17,21,22)]
head(exons.vM20.annotation_gene_chr,11)
head(exons.vM20.annotation_gene_chr)
#Get  exons per gene only one time
exons.vM20.annotation_gene_chr["geneexon_id"] <- data.frame(paste0(exons.vM20.annotation_gene_chr$seqnames, "%",
                                                                   exons.vM20.annotation_gene_chr$start, "%",
                                                                   exons.vM20.annotation_gene_chr$end,"%",
                                                                   exons.vM20.annotation_gene_chr$strand, "%",
                                                                   exons.vM20.annotation_gene_chr$gene_name, "%",
                                                                   exons.vM20.annotation_gene_chr$exon_id))
head(exons.vM20.annotation_gene_chr)
write.table(exons.vM20.annotation_gene_chr, "exons.vM20.annotation_gene_chr.txt", sep = "\t", quote = F, append = F, row.names = F)
sort -k12,12 -u exons.vM20.annotation_gene_chr.txt > exons.vM20.annotation_gene_chro.txt
exons.vM20.annotation_gene_chro <- read.table("exons.vM20.annotation_gene_chro.txt",  header = F, stringsAsFactors = F)
head(exons.vM20.annotation_gene_chro)
dim(exons.vM20.annotation_gene_chro)
colnames(exons.vM20.annotation_gene_chro) <- colnames(exons.vM20.annotation_gene_chr)
head(exons.vM20.annotation_gene_chro)
dim(exons.vM20.annotation_gene_chro)
countexons.vM20.annotation_exon <- count(exons.vM20.annotation_gene_chro, "gene_name")
head(countexons.vM20.annotation_exon)
countexons.vM20.annotation_exon <- countexons.vM20.annotation_exon[order(-countexons.vM20.annotation_exon$freq),]
head(countexons.vM20.annotation_exon)
dim(countexons.vM20.annotation_exon)

countexons.vM20.annotation_gene <- count(exons.vM20.annotation_gene_chro, "geneexon_id")
head(countexons.vM20.annotation_gene)
countexons.vM20.annotation_gene <- countexons.vM20.annotation_gene[order(-countexons.vM20.annotation_gene$freq),]
head(countexons.vM20.annotation_gene)
dim(countexons.vM20.annotation_gene)
library(splitstackshape)
countexons.vM20.annotation_gene_resep <- cSplit(countexons.vM20.annotation_gene, "geneexon_id", "%")
head(countexons.vM20.annotation_gene_resep,1)
countexons.vM20.annotation_gene_resep <- countexons.vM20.annotation_gene_resep[,c(2:6,7,1)]
dim(countexons.vM20.annotation_gene_resep)
head(countexons.vM20.annotation_gene_resep)
countexons.vM20.annotation_gene_resep <- countexons.vM20.annotation_gene_resep[order(countexons.vM20.annotation_gene_resep$geneexon_id_1),]
write.table(countexons.vM20.annotation_gene_resep, "countexons.vM20.annotation_gene_resep.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = F)
grep chr countexons.vM20.annotation_gene_resep.txt > countexons.vM20.annotation_gene_resep_chr.txt

#Count SNPs
bedtools intersect -wa -wb -a countexons.vM20.annotation_gene_resep_chr.txt -b /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-GRCm38-M20-overlapped/jf1v2_Snp+chr.GRCm38.mm10.bed > exons.vM20.annotation_gene_jf1v2.txt
exons.vM20.annotation_gene_jf1v2 <- read.table("exons.vM20.annotation_gene_jf1v2.txt", stringsAsFactors = F)
head(exons.vM20.annotation_gene_jf1v2)
colnames(exons.vM20.annotation_gene_jf1v2) <- c("chr", "start", "end", "strand", "GeneID", "ensID","extra", "chrsnp", "startsnp", "endsnp", "SNP")
countexons.vM20.annotation_gene_jf1v2 <- count(exons.vM20.annotation_gene_jf1v2, "GeneID")
head(countexons.vM20.annotation_gene_jf1v2)
dim(countexons.vM20.annotation_gene_jf1v2)
countexons.vM20.annotation_gene_jf1v2 <- countexons.vM20.annotation_gene_jf1v2[order(-countexons.vM20.annotation_gene_jf1v2$freq),]
head(countexons.vM20.annotation_gene_jf1v2)
colnames(countexons.vM20.annotation_gene_jf1v2) <- c("Genes", "Freq")
write.table(countexons.vM20.annotation_gene_jf1v2, "countexons.vM20.annotation_gene_jf1v2.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = T)
basilia_genes_for_dotplot
countexons.vM20.annotation_gene_jf1v2_npcsig <- merge(countexons.vM20.annotation_gene_jf1v2, basilia_genes_for_dotplot, by="Genes", all.y=TRUE)
write.table(countexons.vM20.annotation_gene_jf1v2_npcsig, "countexons.vM20.annotation_gene_jf1v2_npcsig.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = T)

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  End Of Analysis   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#
