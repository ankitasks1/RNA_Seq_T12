#Note: it should be noted that gtf used for bulk were cleaned for scaffolds and then used for allele-specific 
#..as the chromosome name need to be changed from  chrN to N, scaffolds doesnot have chr only alphabets
# Therefore features in bulk are 55604 and allele sp has 54750
#Strandedness info provided by IGA: First strand = Reversely stranded
#Kit: Illumina TruSeq Stranded Total RNA
#feature count option, s = 2
#No need to sort by name as featurecounts will do it by its own ref1) https://www.biostars.org/p/297865/ ref2)https://bioconductor.org/packages/release/bioc/vignettes/Rsubread/inst/doc/ ref3) https://manpages.debian.org/testing/subread/featureCounts.1.en.htmlSubreadUsersGuide.pdf ref4)https://www.biostars.org/p/245826/ ref5)http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf ref6)https://academic.oup.com/bioinformatics/article/30/7/923/232889


#/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -p -s 2 -a ~/ref_av/gencodes/gencode_M20/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf -T 12 -o JB1_time12_star-featureCounts_mm10.txt JB1_WT_Rep1_Aligned.sortedByCoord.out.bam JB1_WT_Rep2_Aligned.sortedByCoord.out.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bam


#Igf2
#/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -p -s 2 -a ~/ref_av/gencodes/gencode_M20/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotationlike-.gtf -T 12 -o JB1_time12_star-featureCounts_mm10_Igf2.txt JB1_WT_Rep1_Aligned.sortedByCoord.out.bam JB1_WT_Rep2_Aligned.sortedByCoord.out.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bam



#SortedbyReadname (Give same results as above commands because featurecounts sort automatically by readname)

#/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -p -s 2 -a ~/ref_av/gencodes/gencode_M20/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf -T 12 -o JB1_time12_star-featureCounts_mm10.sortbyreadname.txt JB1_WT_Rep1_Aligned.sortedByReadname.out.bam JB1_WT_Rep2_Aligned.sortedByReadname.out.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.bam


#Igf2
#/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -p -s 2 -a ~/ref_av/gencodes/gencode_M20/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotationlike-.gtf -T 12 -o JB1_time12_star-featureCounts_mm10_Igf2.sortbyreadname.txt JB1_WT_Rep1_Aligned.sortedByReadname.out.bam JB1_WT_Rep2_Aligned.sortedByReadname.out.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.bam

##################################################################################### DESeq2 ###############################################################################################################
library(DESeq2)
library(EDASeq)
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2")
countdata <- read.table("JB1_time12_star-featureCounts_mm10_Igf2.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
countdata <- countdata[ ,6:ncol(countdata)]
head(countdata)
# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))
dim(countdata)
head(countdata)
countData = as.matrix(countdata)
head(countData)
colnames(countData) <- c("JB1_WT_Rep1", "JB1_WT_Rep2", "JB1_ZFP57_KO_Rep1", "JB1_ZFP57_KO_Rep2")
coldata <- read.table("coldata.txt" , header = TRUE, stringsAsFactors = FALSE)
rownames(coldata) <- colnames(countData)
head(coldata)
coldata <- coldata[,c("condition","replicate")]
coldata$condition <- factor(coldata$condition)
coldata$replicate <- factor(coldata$replicate)
all(rownames(coldata) == colnames(countData)) #should print TRUE
dds <- DESeqDataSetFromMatrix(countData =countData, colData = coldata, design = ~ condition)
dds
featureData <- data.frame(gene=rownames(countData))
keep <- rowSums(counts(dds)) >= 9
dds <- dds[keep,]
ddscounts <- counts(dds, normalized=FALSE)
head(ddscounts)
dim(ddscounts)
colSums(ddscounts)
#View filtered count matrix: View(counts(dds))
#Normalization is the part of DESeq command: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#Normalized separately:
ddsNorm <- estimateSizeFactors(dds)
sizeFactors(ddsNorm)
#Export normalized counts
ddsNormcounts <- counts(ddsNorm, normalized=TRUE)
head(ddsNormcounts)
dim(ddsNormcounts)
write.table(ddsNormcounts, "ddsNormcounts.txt", sep="\t", quote=F, col.names=NA)
plotPCA(as.matrix(ddsNormcounts), labels=F, col =  c("darkgreen","darkgreen","darkred","darkred"))
#save as PCA_ddsNormcounts.svg

#DEG
#Normalization is the part of DESeq command: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
dds <- DESeq(dds)
#See size factors estimated by DESeq function: sizeFactors(dds)
res <- results(dds, contrast=c("condition","ZFP57KO","WT"))
res
summary(res)
res_sorted = res[order(rownames(res)),]
head(res_sorted)
sum(res$padj < 0.05, na.rm=TRUE)
res_sorted = res[order(rownames(res)),]
head(res_sorted)
summary(res_sorted)
write.table(res_sorted, "deseq2_results_res_sorted.txt", sep="\t", quote = FALSE, append = FALSE)
sum(res_sorted$padj < 0.05, na.rm=TRUE)
res_sorted$threshold <- as.logical(res_sorted$padj < 0.05)
deseq2_results_res0.05 <- res_sorted[which(res_sorted$threshold == TRUE),]
head(deseq2_results_res0.05)
dim(deseq2_results_res0.05)
summary(deseq2_results_res0.05)
write.table(deseq2_results_res0.05, "deseq2_results_res0.05_sorted.txt", sep="\t", quote = FALSE, append = FALSE)

#Tag chromosomes
#head(res_sorted)
res_Deseq2_features <- data.frame(rownames(res_sorted))
colnames(res_Deseq2_features) <- "id"
head(res_Deseq2_features)
#Chromosome positions
chr.pos =  read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt",header=FALSE)
head(chr.pos)
colnames(chr.pos) <- c("id", "Genes", "chr", "start", "end")
head(chr.pos)
chr.pos.1 = merge(res_Deseq2_features, chr.pos, by="id", all.x=TRUE)
head(chr.pos.1)
dim(chr.pos.1)
head(res_sorted)
res_sorted_gene <- data.frame(cbind(res_sorted,chr.pos.1))#just reversed the cbinding to check if row merge exact match
head(res_sorted_gene)
tail(res_sorted_gene)
dim(res_sorted_gene)
#res_sorted_gene_rearranged <- res_sorted_gene[,c(3:5, 1:2, 6:12)]
res_sorted_gene_rearranged <- res_sorted_gene[,c(10:12, 8:9, 1:7)]
res_sorted_gene_rearranged <- res_sorted_gene_rearranged[order(res_sorted_gene_rearranged$chr, res_sorted_gene_rearranged$start),]
head(res_sorted_gene_rearranged)
dim(res_sorted_gene_rearranged)
write.table(res_sorted_gene_rearranged, "res_sorted_gene_rearranged.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F, col.names = F)

##### FOR NOW: Do not filter deletion regions ##############
#Filter non-windows/biased regions : For detailed description see Oct2019 in documents folder 
#Windows identified by allele-specific
#cat NonwINDOWS_matRatioWTsort.txt NonwINDOWS_matRatioZFP57KOsort.txt > combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort.txt
#awk '{if($5 >= 0) print}' combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort.txt > combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort_coverage.txt (38)
#awk '{if($5 < 0) print}' combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort.txt > combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort_NOcoverage.txt (136)
#awk '{print $1"%"$2"%"$3}' combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort_coverage.txt | sort -k1,1 -u | awk -F'%' '{print $1"\t"$2"\t"$3}' | sort -k1,1 -k2,2n > combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort_coverage_toremove.txt (27)
#Check the loss: bedtools intersect -wa -wb -a res_sorted_gene_rearranged.txt -b combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort_coverage_toremove.txt | head
#bedtools intersect -wa -wb -a res_sorted_gene_rearranged.txt -b combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort_coverage_toremove.txt -v > res_sorted_gene_rearranged_filtered.txt

#Tag chromosomes for Deregulated genes
#head(deseq2_results_res0.05)
res_DE_Deseq2_features <- data.frame(rownames(deseq2_results_res0.05))
colnames(res_DE_Deseq2_features) <- "id"
head(res_DE_Deseq2_features)
#Chromosome positions
chr.pos.de = merge(res_DE_Deseq2_features, chr.pos, by="id", all.x=TRUE)
head(chr.pos.de)
dim(chr.pos.de)
head(deseq2_results_res0.05)
deseq2_results_res0.05_gene <- data.frame(cbind(deseq2_results_res0.05,chr.pos.de))
head(deseq2_results_res0.05_gene)
tail(deseq2_results_res0.05_gene)
dim(deseq2_results_res0.05_gene)
deseq2_results_res0.05_gene_rearranged <- deseq2_results_res0.05_gene[,c(10:12, 8:9, 1:7)]
deseq2_results_res0.05_gene_rearranged <- deseq2_results_res0.05_gene_rearranged[order(deseq2_results_res0.05_gene_rearranged$chr, deseq2_results_res0.05_gene_rearranged$start),]
head(deseq2_results_res0.05_gene_rearranged)
dim(deseq2_results_res0.05_gene_rearranged)
write.table(deseq2_results_res0.05_gene_rearranged, "deseq2_results_res0.05_gene_rearranged.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F, col.names = F)


#Volcano Plot extra
library("dplyr")
library("ggplot2")
res_sorted_gene_rearranged_filtered <- read.table("res_sorted_gene_rearranged.txt", header = F) #Do not filter windows for now
head(res_sorted_gene_rearranged_filtered)
colnames(res_sorted_gene_rearranged_filtered) <- c("chr", "start", "end", "ensID", "Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "threshold")
head(res_sorted_gene_rearranged_filtered)
dim(res_sorted_gene_rearranged_filtered)
res_sorted_gene_WT_KO <- res_sorted_gene_rearranged_filtered[,c(5,7,11)]
head(res_sorted_gene_WT_KO)
results = mutate(res_sorted_gene_WT_KO, sig=ifelse(res_sorted_gene_WT_KO$padj<0.05, "Padj<0.05", "Not_Sig"))
head(results)
write.table(results, "results_res_sorted_gene_WT_KO.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)

# Make a basic volcano plot #https://www.gettinggeneticsdone.com/2014/05/r-volcano-plots-to-visualize-rnaseq-microarray.html
with(results, plot(log2FoldChange, -log10(padj), pch=21, cex = 0.3, main="Volcano plot", xlim=c(-15,15), col="black", ylim=c(0,100)))
#Coloroffsets
#Upregulated
with(subset(results, padj<.05 & log2FoldChange < 1), points(log2FoldChange, -log10(padj), pch=21, cex = 0.4, xlim=c(-15,15), ylim=c(0,100),lwd = 0.4, col="black",bg="grey"))
#Downregulated
with(subset(results, padj<.05 & log2FoldChange > -1), points(log2FoldChange, -log10(padj), pch=21, cex = 0.4, xlim=c(-15,15), ylim=c(0,100), lwd = 0.4,col="black", bg="grey"))

#Upregulated
with(subset(results, padj<.05 & log2FoldChange > 1), points(log2FoldChange, -log10(padj), pch=21, cex = 0.4, xlim=c(-15,15), ylim=c(0,100),lwd = 0.4, col="darkgreen",bg="#66CDAA"))
#Downregulated
with(subset(results, padj<.05 & log2FoldChange< -1), points(log2FoldChange, -log10(padj), pch=21, cex = 0.4, xlim=c(-15,15), ylim=c(0,100), lwd = 0.4,col="darkred", bg="#F08080"))
#Imprinted genes p.adj <0.05
#fgrep -f imprinted_gene_name.txt results_res_sorted_gene_WT_KO.txt -w | grep Padj > results_imprinted.txt
results_imprinted <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/results_imprinted.txt",header=FALSE)
colnames(results_imprinted) <- c("Genes",	"log2FoldChange",	"padj",	"sig")
head(results_imprinted)
dim(results_imprinted)
with(subset(results_imprinted, padj<.05 & log2FoldChange > 0), points(log2FoldChange, -log10(padj), pch=21, cex = 0.6, xlim=c(-15,15), ylim=c(0,100), lwd = 1,col="orange", bg="yellow"))
with(subset(results_imprinted, padj<.05 & log2FoldChange < 0), points(log2FoldChange, -log10(padj), pch=21, cex = 0.6, xlim=c(-15,15), ylim=c(0,100), lwd = 1,col="orange", bg="yellow"))

#ICR genes p.adj <0.05
#fgrep -f icr.genes results_res_sorted_gene_WT_KO.txt -w | grep Padj > results_icr.txt
results_icr <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/results_icr.txt",header=FALSE)
colnames(results_icr) <- c("Genes",	"log2FoldChange",	"padj",	"sig")
head(results_icr)
dim(results_icr)
with(subset(results_icr, padj<.05 & log2FoldChange > 0), points(log2FoldChange, -log10(padj), pch=21, cex = 0.7, xlim=c(-15,15), ylim=c(0,100), lwd = 1,col="blue", bg="#87CEFA"))
with(subset(results_icr, padj<.05 & log2FoldChange < 0), points(log2FoldChange, -log10(padj), pch=21, cex = 0.7, xlim=c(-15,15), ylim=c(0,100), lwd = 1,col="blue", bg="#87CEFA"))

#text(results_icr$log2FoldChange, -log10(results_icr$padj), labels=results_icr$Genes, cex= 0.4, pos=3, col = "black")
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
#with(subset(results, padj<.05 ), points(log2FoldChange, -log10(padj), pch=20, cex = 0.5, xlim=c(-20,20), ylim=c(0,100), col="black"))
#with(subset(results, abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, cex = 0.5, xlim=c(-20,20), ylim=c(0,100), col="orange"))
#with(subset(results, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, cex = 0.5, xlim=c(-20,20), ylim=c(0,100), col="green"))

#Gene lenth supplied in RPKM is the sum of length of exons,
#as in this example,: sum(c(259, 43, 142, 158, 129, 130, 154, 210, 2037)) # see Ensembl for gene lengths,
#This gene length is already given by featurecounts , so use it.
countforfpkm <- read.table("JB1_time12_star-featureCounts_mm10_Igf2.txt", header=TRUE, row.names=1)
countforfpkm <- countforfpkm[ ,5:ncol(countforfpkm)]
head(countforfpkm)
countforfpkm <- countforfpkm[order(rownames(countforfpkm)),]
head(countforfpkm)
countforfpkm_filt <- countforfpkm[which(rowSums(countforfpkm[,2:5]) >=0),] #Dont filter for now use all the genes
dim(countforfpkm_filt)
head(countforfpkm_filt)

gene.length <- countforfpkm_filt$Length
countforfpkm_filt1 <- countforfpkm_filt[,2:5]
head(countforfpkm_filt1)
colnames(countforfpkm_filt1) <-  c("JB1_WT_Rep1","JB1_WT_Rep2","JB1_ZFP57_KO_Rep1","JB1_ZFP57_KO_Rep2")
head(countforfpkm_filt1)

#Use edgeR for FPKM calculation, the method example is (readcount * 10^9)/(Depth * readlength * normfactor) eg. (1701 * 10^9) / (7661385 * 3262 * 1.0577571) = 64.34682
#edgeR rpkm()
library(edgeR)
dgeforfpkm <- DGEList(counts=countforfpkm_filt1,genes=data.frame(Length=gene.length), group = c("WT", "WT", "ZFP57KO", "ZFP57KO"))
dgeforfpkm <- calcNormFactors(dgeforfpkm)
dgeFPKM <- rpkm(dgeforfpkm, dgeforfpkm$genes$Length)
head(dgeFPKM)
write.table(dgeFPKM, "dgeFPKM.txt", sep="\t", quote = FALSE, append = FALSE)

#Extract  markers:
#fgrep -f selected.markers /home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20_rearranged.txt -w | awk '{print $4}' | sort -k1,1 > selected.markers.ens
#fgrep -f selected.markers /home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20_rearranged.txt -w | sort -k4,4 > selected.markers.ensgene
#fgrep -f selected.markers.ens dgeFPKM.txt -w | sort -k1,1 > dgeFPKM.selectedmarkers.txt
#paste selected.markers.ensgene dgeFPKM.selectedmarkers.txt > marker_time12_heatmap.txt

marker_time12_heatmap <- read.table("marker_time12_heatmap.txt", header = F)
head(marker_time12_heatmap)
marker_time12_heatmap1 <- marker_time12_heatmap[,c(5,7:10)]
head(marker_time12_heatmap1)
colnames(marker_time12_heatmap1) <- c("Gene", "WT1", "WT2", "ZFP57KO1", "ZFP57KO2")
head(marker_time12_heatmap1)
write.table(marker_time12_heatmap1,"marker_time12_heatmap1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)

#or directly import T0 data and combine with T12
#Heatmap
library(ggplot2)
library(gplots)
dgeFPKM <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/dgeFPKM.txt")
#T12 FPKM
dgeFPKM_T12 <- dgeFPKM
dim(dgeFPKM_T12)
head(dgeFPKM_T12)
colnames(dgeFPKM_T12) <- c("WT1.t12", "WT2.t12", "ZFP57KO1.t12","ZFP57KO2.t12")
dgeFPKM_T0 <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2018/time0/2015_05_28_JB1_WT_Zfp57--_clone4/DATI_GREZZI/Bulk_mm10/featurecount/Igf2/dgeFPKM.txt", header = T)
head(dgeFPKM_T0)
dim(dgeFPKM_T0)
colnames(dgeFPKM_T0) <- c("WT1.t0", "WT2.t0", "WT3.t0", "ZFP57KO1.t0","ZFP57KO2.t0","ZFP57KO3.t0")
dgeFPKM_T12_T0 <-  cbind.data.frame(rownames(dgeFPKM_T12), dgeFPKM_T12, rownames(dgeFPKM_T0),dgeFPKM_T0)
head(dgeFPKM_T12_T0)
tail(dgeFPKM_T12_T0)
dim(dgeFPKM_T12_T0)
#average samples by rowmeans, I manually checked this commands
dgeFPKM_T12_T0_avg <- cbind.data.frame(rownames(dgeFPKM_T12_T0),
                                       rowMeans(dgeFPKM_T12_T0[,2:3]),
                                       rowMeans(dgeFPKM_T12_T0[,4:5]),
                                       rowMeans(dgeFPKM_T12_T0[,7:9]),
                                       rowMeans(dgeFPKM_T12_T0[,10:12]))
colnames(dgeFPKM_T12_T0_avg) <- c("id", "WT_T12", "ZFP57KO_T12", "WT_T0", "ZFP57KO_T0")
head(dgeFPKM_T12_T0_avg)
tail(dgeFPKM_T12_T0_avg)

#Average only T0 and keep T12 sep
head(dgeFPKM_T12_T0)
dgeFPKM_T12_T0_t0avg <- cbind.data.frame(rownames(dgeFPKM_T12_T0),
                                       (dgeFPKM_T12_T0[,2:3]),
                                       (dgeFPKM_T12_T0[,4:5]),
                                       rowMeans(dgeFPKM_T12_T0[,7:9]),
                                       rowMeans(dgeFPKM_T12_T0[,10:12]))
head(dgeFPKM_T12_T0_t0avg)
colnames(dgeFPKM_T12_T0_t0avg) <- c("id", "WT1_T12", "WT2_T12", "ZFP57KO1_T12", "ZFP57KO2_T12", "WT_T0", "ZFP57KO_T0")
head(dgeFPKM_T12_T0_t0avg)
tail(dgeFPKM_T12_T0_t0avg)
#Assign chromsomes and gene names
head(chr.pos)
dgeFPKM_T12_T0_t0avg.chr = merge(dgeFPKM_T12_T0_t0avg, chr.pos, by="id", all.x=TRUE)
head(dgeFPKM_T12_T0_t0avg.chr)
dgeFPKM_T12_T0_t0avg.chr.gene <- dgeFPKM_T12_T0_t0avg.chr[,c(8,6,7,2:5)]
head(dgeFPKM_T12_T0_t0avg.chr.gene)
dim(dgeFPKM_T12_T0_t0avg.chr.gene)
write.table(dgeFPKM_T12_T0_t0avg.chr.gene,"dgeFPKM_T12_T0_t0avg.chr.gene.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)

#selected_markers_T0_T12 <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/heatmap-time0-time12_markers.txt", header = T)
#rownames(selected_markers_T0_T12)
selected.markers <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/selected.markers", header = F)
colnames(selected.markers) <- "Genes"
selected_markers_T0_T12_1 = merge(selected.markers, dgeFPKM_T12_T0_t0avg.chr.gene, by="Genes", all.x=TRUE, sort = FALSE)
head(selected_markers_T0_T12_1)
dim(selected_markers_T0_T12_1)
selected_markers_T0_T12_1 <- selected_markers_T0_T12_1[1:25,]
rownames(selected_markers_T0_T12_1)=selected_markers_T0_T12_1[,1]
rownames(selected_markers_T0_T12_1)
colnames(selected_markers_T0_T12_1)
head(selected_markers_T0_T12_1)
selected_markers_T0_T12_1 = selected_markers_T0_T12_1[,-1]
selected_markers_T0_T12_1 = as.matrix(selected_markers_T0_T12_1)
head(selected_markers_T0_T12_1)
dim(selected_markers_T0_T12_1)
write.table(selected_markers_T0_T12_1,"/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/selected_markers_T0_T12_1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)
selected_markers_T0_T12_1 <- data.frame(selected_markers_T0_T12_1)

#Beta-actin based normalization Actb	328.893776284484, 612.889483712133, 673.968354045603, 1348.52877698723, 1396.11296026648, 1501.60527280486
head(selected_markers_T0_T12_1)
grep Actb dgeFPKM_T12_T0_t0avg.chr.gene.txt
#Actb	328.893776284484	612.889483712133	673.968354045603	1348.52877698723	1396.11296026648	1501.60527280486

Actb <- c(328.893776284484, 612.889483712133, 673.968354045603, 1348.52877698723, 1396.11296026648, 1501.60527280486)
Actb <- as.matrix(Actb)
selected_markers_T0_T12_1 <- as.matrix(selected_markers_T0_T12_1)
selected_markers_T0_T12_1.norm <- cbind.data.frame((selected_markers_T0_T12_1[,1] / Actb[1]),
                                                   (selected_markers_T0_T12_1[,2] / Actb[2]),
                                                   (selected_markers_T0_T12_1[,3] / Actb[3]),
                                                   (selected_markers_T0_T12_1[,4] / Actb[4]),
                                                   (selected_markers_T0_T12_1[,5] / Actb[5]),
                                                   (selected_markers_T0_T12_1[,6] / Actb[6]))
head(selected_markers_T0_T12_1.norm,1)
colnames(selected_markers_T0_T12_1.norm) <- c("WT_T0", "ZFP57KO_T0", "WT1_T12", "WT2_T12", "ZFP57KO1_T12", "ZFP57KO2_T12")
head(selected_markers_T0_T12_1.norm)
write.table(selected_markers_T0_T12_1.norm,"/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/selected_markers_T0_T12_1.norm.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)

selected_markers_T0_T12_1.norm <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/selected_markers_T0_T12_1.norm.txt", header = T)
head(selected_markers_T0_T12_1.norm)
dim(selected_markers_T0_T12_1.norm)

selected_markers_T0_T12_1.norm.re <- selected_markers_T0_T12_1.norm[c(1:12,15,17:25),]
head(selected_markers_T0_T12_1.norm.re)
# load package
library(pheatmap)
library(RColorBrewer)
breaksList = seq(0, 0.1, by = 0.001)
my_sample_row1 <- data.frame(Markers= c("ESC","ESC","ESC","ESC","ESC","Neural","Neural","Neural","Neural","Neural","Neural","Neural","Cortical","Cortical","Endo_Meso","Endo_Meso","Endo_Meso","Endo_Meso","Endo_Meso","Endo_Meso","Endo_Meso","Knockout"))
rownames(my_sample_row1) <- rownames(selected_markers_T0_T12_1.norm.re)
my_colour1 = list(Markers = c(ESC = "#e6e600", Neural = "#00e6ac", Cortical ="#E0BBE4", Endo_Meso = "#ff80b3", Knockout = "#B3A580"))
#pheatmap(selected_markers_T0_T12_1.ratio,treeheight_row = 0, cluster_cols=F, cluster_rows=T, treeheight_col = 0, gaps_col =NULL, gaps_row = NULL, border_color = "black", breaks = breaksList,color= colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)))
#Remove Sox1
pheatmap(selected_markers_T0_T12_1.norm.re[c(1:5,7:22),],
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(length(breaksList)),
         breaks = breaksList,
         annotation_colors = my_colour1,
         annotation_row = my_sample_row1,
         cluster_cols=F, 
         cluster_rows=F,
         fontsize = 12,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cutree_cols = 1)

#save as selected_markers_T0_T12_1.norm.re.svg

colfunc <- colorRampPalette(c("#2166AC","white","#B2182B"))
selected_markers_T0_T12_1.norm.re <- as.matrix(selected_markers_T0_T12_1.norm.re)

heatmap.2(selected_markers_T0_T12_1.norm.re[c(1:5,7:22),], 
          Colv = "NA", 
          Rowv = "NA",
          trace = "none", 
          col = colfunc , 
          lmat=rbind( c(5, 4, 2), c(6,1,3)), 
          lwid=c(1, 4,2),
          lhei = c(1,7), 
          keysize=1.2, 
          key.par = list(cex=0.6), 
          density.info=c("none"), 
          dendrogram="none", 
          scale = "none", 
          sepwidth=c(0.001, 0.001), 
          cexRow=2, font=3, 
          cexCol = 1, 
          margins =c(5,8), 
          breaks = seq(0,0.2, length.out = 100))

#Scatter:Log normalized counts
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/")
ddsNormcounts <- read.table("ddsNormcounts.txt" , header = TRUE, stringsAsFactors = FALSE)
head(ddsNormcounts)
dim(ddsNormcounts)
colnames(ddsNormcounts) <- c("JB1_WT_Rep1","JB1_WT_Rep2","JB1_ZFP57_KO_Rep1","JB1_ZFP57_KO_Rep2")
ddsNormcounts <- data.frame(ddsNormcounts[order(rownames(ddsNormcounts)),])
head(ddsNormcounts)
ddsnorm_DE_Deseq2_features <- data.frame(rownames(ddsNormcounts))
colnames(ddsnorm_DE_Deseq2_features) <- "id"
head(ddsnorm_DE_Deseq2_features)
#Chromosome positions
chr.pos =  read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt",header=FALSE)
head(chr.pos)
colnames(chr.pos) <- c("id", "Genes", "chr", "start", "end")
head(chr.pos)
chr.pos.1 = merge(ddsnorm_DE_Deseq2_features, chr.pos, by="id", all.x=TRUE)
head(chr.pos.1)
tail(chr.pos.1)
dim(chr.pos.1)
head(ddsNormcounts)
tail(ddsNormcounts)
ddsNormcounts_gene <- data.frame(cbind(ddsNormcounts,chr.pos.1))
head(ddsNormcounts_gene)
tail(ddsNormcounts_gene)
write.table(ddsNormcounts_gene, "ddsNormcounts_gene.txt", sep="\t", quote = FALSE, append = FALSE, row.names = TRUE)
ddsNormcountsavg <- data.frame(cbind(ddsNormcounts_gene[,c(5:9)],(ddsNormcounts_gene$JB1_WT_Rep1+ddsNormcounts_gene$JB1_WT_Rep2)/2, (ddsNormcounts_gene$JB1_ZFP57_KO_Rep1 + ddsNormcounts_gene$JB1_ZFP57_KO_Rep2)/2))
head(ddsNormcountsavg)
#rownames(ddsNormcountsavg) <- rownames(ddsNormcounts)
dim(ddsNormcountsavg)
colnames(ddsNormcountsavg) <- c("id","Genes","chr","start","end","JB1_WT","JB1_ZFP57KO")
rownames(ddsNormcountsavg)
head(ddsNormcountsavg)
tail(ddsNormcountsavg)
write.table(ddsNormcountsavg, "ddsNormcountsavg_gene_WT_KO.pre.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)

ddsNormcountsavgLog <- data.frame(cbind(ddsNormcountsavg[,c(3:5,2,1)], log2(ddsNormcountsavg[,c(6:7)]))) #Take log2 of DESeq2 normalized count
head(ddsNormcountsavgLog)
dim(ddsNormcountsavgLog)
write.table(ddsNormcountsavgLog, "ddsNormcountsavgLog_gene_WT_KO.pre.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)
#Filter non-windows/biased regions : For detailed description see Oct2019 in documents folder 
#grep ENS ddsNormcountsavgLog_gene_WT_KO.pre.txt | sort -k1,1 -k2,2n > ddsNormcountsavgLog_gene_WT_KO.sorted.txt
#ddsNormcountsavgLog_gene_WT_KO.sorted.txt
###Don't run this command, Do not filter for now #####bedtools intersect -wa -wb -a ddsNormcountsavgLog_gene_WT_KO.sorted.txt -b combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort_coverage_toremove.txt -v > ddsNormcountsavgLog_gene_WT_KO.sorted1.txt 
###Don't run this command, Do not filter for now ##### cp ddsNormcountsavgLog_gene_WT_KO.sorted1.txt all_ddsNormcountsavgLog_gene_WT_KO.sorted.txt

cp ddsNormcountsavgLog_gene_WT_KO.sorted.txt all_ddsNormcountsavgLog_gene_WT_KO.sorted.txt
awk '{print $1}' deseq2_results_res0.05_sorted.txt | grep ENS > deseq2_results_res0.05_sorted.ensid

#188 imprinted_gene_name.txt
fgrep -f imprinted_gene_name.txt all_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w > imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.txt 
fgrep -f imprinted_gene_name.txt all_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w -v > Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.txt
#Deregulated
fgrep -f deseq2_results_res0.05_sorted.ensid all_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w > all_ddsNormcountsavgLog_gene_WT_KO.sorted.deregulated.txt
fgrep -f deseq2_results_res0.05_sorted.ensid imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w > imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.deregulated.txt
fgrep -f deseq2_results_res0.05_sorted.ensid Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w > Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.deregulated.txt

#Non deregulated
fgrep -f deseq2_results_res0.05_sorted.ensid all_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w -v > all_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.txt
fgrep -f deseq2_results_res0.05_sorted.ensid imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w -v > imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.txt
fgrep -f deseq2_results_res0.05_sorted.ensid Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w -v > Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.txt

#For more commands to generate scatter file follow October2019_1-15.odt
#...................
#..................
#Peak calling:
#ZFP57: macs14 -t SRR2927818_rmdup.sorted1.bed  -c ./../Input/SRR2927817_rmdup.sorted1.bed -g mm -n peak_SRR2927817_818 -B -S --call-subpeaks
#KAP1: macs14 -t SRR2927819_rmdup.sorted1.bed  -c ./../Input/SRR2927817_rmdup.sorted1.bed -g mm -n peak_SRR2927817_819 -B -S --call-subpeaks

#Overlap JB1_Zfp57 (544) and JB1_Kap1 (6986) peaks and extract only JB1_Zfp57 peaks (457)
#bedtools intersect -wa -wb -a Zfp57_GSM1931786_mm10_peaks.bed -b Kap1_GSM1931787_mm10_peaks.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > Zfp57_overlapped_Kap1_mm10_peaks.bed

#This will generate scatter_bulk_mm10_t12_column.txt file
#Replace color name with numbers and remove -inf containing genes and save as  scatter_bulk_mm10_t12_col_num.txt:I dis it with R so not manually
#awk '{print $1}' scatter_bulk_mm10_t12_col_num.txt > scatter_bulk_mm10_t12_col_num.gene
#library(ggplot2)
#library(ggrepel)
#scatter <- read.table("scatter_bulk_mm10_t12_col_num.txt",header=T)
#rownames(scatter)
#rownames(scatter)
#colnames(scatter)
#head(scatter)
#scatter1 = scatter
#head(scatter1)
#dim(scatter1)
#colnames(scatter1) <- c("Gene", "WT","ZFP57KO","Color")
#head(scatter1)
#with(scatter1, plot(WT, ZFP57KO, pch=21, cex = 0.3, main="Scatter Plot", xlim=c(0,15), ylim=c(0,15),lwd = 0.3,col="#525254", bg="grey",ylab="JB1 Zfp57-/- log2(normalized read counts)", xlab="JB1 WT log2(normalized read counts)", bty = 'n'))
#with(subset(scatter1, Color==3), points(WT, ZFP57KO, pch=21, cex = 0.5, xlim=c(0,15), ylim=c(0,15), lwd = 0.5,col="#bf2828", bg="#fcdcdc"),ylab="", xlab="", bty = 'n') #reddish
#with(subset(scatter1, Color==4), points(WT, ZFP57KO, pch=21, cex = 0.5, xlim=c(0,15), ylim=c(0,15), lwd = 0.5,col="#ffee29", bg="#fffbcc"),ylab="", xlab="", bty = 'n') #yellowish
#with(subset(scatter1, Color==1), points(WT, ZFP57KO, pch=21, cex = 0.8, xlim=c(0,15), ylim=c(0,15), lwd = 1,col="#159632", bg="#c9ffd5"),ylab="", xlab="", bty = 'n') #greenish
#with(subset(scatter1, Color==2), points(WT, ZFP57KO, pch=21, cex = 0.8, xlim=c(0,15), ylim=c(0,15), lwd = 1,col="#1b38b5", bg="#d4dcff"),ylab="", xlab="", bty = 'n') #bluish
#with(subset(scatter1, Color==9), points(WT, ZFP57KO, pch=21, cex = 1.5, xlim=c(0,15), ylim=c(0,15), lwd = 1,col="black", bg="#A0522D"),ylab="JB1 Zfp57-/- log2(normalized read counts)", xlab="JB1 WT log2(normalized read counts)", bty = 'n') #brownish
#lines(x = c(0,15), y = c(0,15), col = "#464647")
#Save manually as scatter_bulk_mm10_t12_col_num-inf.svg


#Overlap JB1_Zfp57 (544) and JB1_Kap1 (6986) peaks and extract only JB1_Zfp57 peaks (unique 457)
bedtools intersect -wa -wb -a Zfp57_GSM1931786_mm10_peaks.bed -b Kap1_GSM1931787_mm10_peaks.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > Zfp57_overlapped_Kap1_mm10_peaks.bed


awk '{if($3>0) print $0}' deseq2_results_res0.05_sorted.txt | awk '{print $1}' | grep ENS > deseq2_results_res0.05_sorted.upregulated.ensid
awk '{if($3<0) print $0}' deseq2_results_res0.05_sorted.txt | awk '{print $1}' | grep ENS > deseq2_results_res0.05_sorted.downregulated.ensid
grep ENS ddsNormcountsavgLog_gene_WT_KO.pre.txt | sort -k1,1 -k2,2n > all_ddsNormcountsavgLog_gene_WT_KO.sorted.txt

#188 imprinted_gene_name.txt
fgrep -f imprinted_gene_name.txt all_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w > imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.txt #123
fgrep -f imprinted_gene_name.txt all_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w -v > Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.txt #19872

#upregulated
fgrep -f deseq2_results_res0.05_sorted.upregulated.ensid all_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w | sort -k1,1 -k2,2n > all_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.txt #2243
fgrep -f deseq2_results_res0.05_sorted.upregulated.ensid imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w | sort -k1,1 -k2,2n > imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.txt #23
fgrep -f deseq2_results_res0.05_sorted.upregulated.ensid Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w | sort -k1,1 -k2,2n > Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.txt #2220

#downregulated
fgrep -f deseq2_results_res0.05_sorted.downregulated.ensid all_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w | sort -k1,1 -k2,2n > all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.txt #2669
fgrep -f deseq2_results_res0.05_sorted.downregulated.ensid imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w | sort -k1,1 -k2,2n > imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.txt #27
fgrep -f deseq2_results_res0.05_sorted.downregulated.ensid Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w | sort -k1,1 -k2,2n > Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.txt #2642


#Non-deregulated
fgrep -f deseq2_results_res0.05_sorted.ensid all_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w -v > all_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.txt #15083
fgrep -f deseq2_results_res0.05_sorted.ensid imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w -v > imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.txt #73
fgrep -f deseq2_results_res0.05_sorted.ensid Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.txt -w -v > Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.txt #15010

#Counts
19995 all_ddsNormcountsavgLog_gene_WT_KO.sorted.txt
123 imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.txt #Out of 188 only 123 overlapped, rest are filtered during low count filtering
19872 Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.txt
#All deregulated  =up-2243 + down-2669 =4912
#Upregulated,2243
2243 all_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.txt
23 imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.txt
2220 Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.txt

#Downregulated,2669
2669 all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.txt
27 imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.txt
2642 Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.txt

#Nonderegulated =19995-4912 =15083
15083 all_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.txt
73 imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.txt
15010 Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.txt

#Possible Effect of Zfp57-/- on gene regulation #Some genes gets repated due to more than one Zfp57 peaks overlapped
##Downregulated
#All (<100kb:2675=188+2487)
bedtools closest -a all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""downregulatedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""downregulatedNotunder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
#Imprinted: 
bedtools closest -a imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""downregulatedImprintedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""downregulatedImprintedNotUnder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
#Nonimprinted: 
bedtools closest -a Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""downregulatedNotImprintedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""downregulatedNotImprintedNotUnder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt

##Upregulated
#All
bedtools closest -a all_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""upregulatedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > all_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a all_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""upregulatedNotunder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > all_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
#Imprinted: 
bedtools closest -a imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""upregulatedImprintedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""upregulatedImprintedNotUnder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
#Nonimprinted: 
bedtools closest -a Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""upregulatedNotImprintedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""upregulatedNotImprintedNotUnder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt

#Nonderegulated
#All
bedtools closest -a all_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""NonDeregulatedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > all_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a all_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""NonDeregulatedNotunder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > all_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
#Imprinted: 
bedtools closest -a imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""NonDeregulatedImprintedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""NonDeregulatedImprintedNotUnder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
#Nonimprinted: 
bedtools closest -a Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""NonDeregulatedNotImprintedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""NonDeregulatedNotImprintedNotUnder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
#count.sup2.sh
Dowregualted (original vs multiple repeated due to intersection with Zfp57peaks)
All (2669 vs 2675) 6
2675 = 188 + 2487
Imprinted
30 =15+15 (27 vs 30) 3
Nonimprinted (2642 vs 2645) 3
2645 =173+2472
Upregulated
All (2243 vs 2246) 3
2246 =160+2086
Imprinted (23 vs 25)2
25=13+12
Nonimprinted (2220 vs 2221)1
2221=147 + 2074
#upper calculation match with down one when run separtely <=100 and >100kb./count.sup.sh
Dowregualted
All
188
2487
Imprinted
15
15
Nonimprinted
173
2472
Upregulated
All
160
2086
Imprinted
13
12
Nonimprinted
147
2074

wc -l all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l all_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l all_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l all_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l all_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt

#save as ./count.sh
188 all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
2487 all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
15 imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
15 imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
173 Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
2472 Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
160 all_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
2086 all_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
13 imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
12 imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
147 Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
2074 Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
1011 all_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
14088 all_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
24 imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
50 imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
987 Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
14038 Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt

#15+15+173+2472+13+12+147+2074+24+50+987+14038 =20020
cat imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt  imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt imprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt  Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt Nonimprinted_ddsNormcountsavgLog_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt | sort -k5,5 > scatterundbulk_mm10_t12.txt

awk '{print $4"\t"$6"\t"$7"\t"$10}'  scatterundbulk_mm10_t12.txt > scatterund_bulk_mm10_t12_column.txt
awk  '{print $4}' scatterund_bulk_mm10_t12_column.txt | sort -k1,1 -u #print all possibilities

downregulatedImprintedNotUnder100Kb 1
downregulatedImprintedUnder100Kb 2
downregulatedNotImprintedNotUnder100Kb 3
downregulatedNotImprintedUnder100Kb 4
NonDeregulatedImprintedNotUnder100Kb 9
NonDeregulatedImprintedUnder100Kb  10
NonDeregulatedNotImprintedNotUnder100Kb 11
NonDeregulatedNotImprintedUnder100Kb  12
upregulatedImprintedNotUnder100Kb  5
upregulatedImprintedUnder100Kb 6
upregulatedNotImprintedNotUnder100Kb  7
upregulatedNotImprintedUnder100Kb 8
Zfp57 13

#Replace effect name with numbers and save as  scatterund_bulk_mm10_t12_col_num.txt #-inf not found
sed 's/downregulatedImprintedNotUnder100Kb/1/g' scatterund_bulk_mm10_t12_column.txt | sed 's/downregulatedImprintedUnder100Kb/2/g' | sed 's/downregulatedNotImprintedNotUnder100Kb/3/g' | sed 's/downregulatedNotImprintedUnder100Kb/4/g' | sed 's/upregulatedImprintedNotUnder100Kb/5/g' | sed 's/upregulatedImprintedUnder100Kb/6/g' | sed 's/upregulatedNotImprintedNotUnder100Kb/7/g' | sed 's/upregulatedNotImprintedUnder100Kb/8/g' | sed 's/NonDeregulatedImprintedNotUnder100Kb/9/g' | sed 's/NonDeregulatedImprintedUnder100Kb/10/g' | sed 's/NonDeregulatedNotImprintedNotUnder100Kb/11/g' | sed 's/NonDeregulatedNotImprintedUnder100Kb/12/g' > scatterund_bulk_mm10_t12_col_num.txt

#Assign ZFP57
#grep Zfp57 scatterund_bulk_mm10_t12_column.txt -w   > Zfp57dot.txt
#Zfp57	8.39246259782179	7.19688329910957	downregulatedNotImprintedNotUnder100Kb
#vim Zfp57dot.txt: downregulatedNotImprintedNotUnder100Kb to 13 and Zfp57 to Zfp57dot
cat scatterund_bulk_mm10_t12_col_num.txt Zfp57dot.txt > scatterund_bulk_mm10_t12_col_num.Zfp57dot.txt
grep Inf2 scatterund_bulk_mm10_t12_col_num.Zfp57dot.txt > Inf2.txt #Keep actual Inf gene symbol
#Remove Inf genes (infinite values)
grep Inf scatterund_bulk_mm10_t12_col_num.Zfp57dot.txt -v > scatterund_bulk_mm10_t12_col_num.Zfp57dot_Inf.txt
#---One imprinted gene cab be plotted due to -Inf, (Kcnk9	7.61445821702791	-Inf	1)
#Re-add Inf gene, gene name is Inf2
cat scatterund_bulk_mm10_t12_col_num.Zfp57dot_Inf.txt Inf2.txt > scatterund_bulk_mm10_t12_col_num.Zfp57dot_Inf2.txt
#Left overs after removing Inf genes = 19454 genes 
library(ggplot2)
library(ggrepel)
scatterund <- read.table("scatterund_bulk_mm10_t12_col_num.Zfp57dot_Inf2.txt",header=F)
head(scatterund)
scatterund_re <- scatterund
rownames(scatterund_re)
colnames(scatterund_re)
head(scatterund_re)
dim(scatterund_re)
write.table(scatterund_re, "scatterund_re.txt", quote = F, append = F)

scatterund1 <- scatterund_re
head(scatterund1)
dim(scatterund1)
colnames(scatterund1) <- c("Gene", "WT","ZFP57KO","Color")
head(scatterund1)
tail(scatterund1)
dim(scatterund1)
scatterund1 <- data.frame(scatterund1)
with(scatterund1, plot(WT, ZFP57KO, pch=21, cex = 0.2, main="Scatter Plot", xlim=c(0,15), ylim=c(0,15),lwd = 0.3,col="#525254", bg="grey",ylab="JB1 Zfp57-/- log2(normalized read counts)", xlab="JB1 WT log2(normalized read counts)", bty = 'n'))
with(subset(scatterund1, Color==7), points(WT, ZFP57KO, pch=21, cex = 0.4, xlim=c(0,15), ylim=c(0,15), lwd = 0.5,col="#9fc5e8ff", bg="#9fc5e8ff"),ylab="", xlab="", bty = 'n') #reddish
with(subset(scatterund1, Color==8), points(WT, ZFP57KO, pch=21, cex = 0.4, xlim=c(0,15), ylim=c(0,15), lwd = 0.5,col="#9fc5e8ff", bg="#9fc5e8ff"),ylab="", xlab="", bty = 'n') #yellowish
with(subset(scatterund1, Color==5), points(WT, ZFP57KO, pch=21, cex = 0.6, xlim=c(0,15), ylim=c(0,15), lwd = 1,col="#0b5394ff", bg="#0b5394ff"),ylab="", xlab="", bty = 'n') #greenish
with(subset(scatterund1, Color==6), points(WT, ZFP57KO, pch=21, cex = 0.6, xlim=c(0,15), ylim=c(0,15), lwd = 1,col="#0b5394ff", bg="#0b5394ff"),ylab="", xlab="", bty = 'n') #bluish
with(subset(scatterund1, Color==3), points(WT, ZFP57KO, pch=21, cex = 0.4, xlim=c(0,15), ylim=c(0,15), lwd = 0.5,col="#f4ccccff", bg="#f4ccccff"),ylab="", xlab="", bty = 'n') #reddish
with(subset(scatterund1, Color==4), points(WT, ZFP57KO, pch=21, cex = 0.4, xlim=c(0,15), ylim=c(0,15), lwd = 0.5,col="#f4ccccff", bg="#f4ccccff"),ylab="", xlab="", bty = 'n') #yellowish
with(subset(scatterund1, Color==1), points(WT, ZFP57KO, pch=21, cex = 0.6, xlim=c(0,15), ylim=c(0,15), lwd = 1,col="#cc0000ff", bg="#cc0000ff"),ylab="", xlab="", bty = 'n') #greenish
with(subset(scatterund1, Color==2), points(WT, ZFP57KO, pch=21, cex = 0.6, xlim=c(0,15), ylim=c(0,15), lwd = 1,col="#cc0000ff", bg="#cc0000ff"),ylab="", xlab="", bty = 'n') #bluish

with(subset(scatterund1, Color==13), points(WT, ZFP57KO, pch=21, cex = 1.2, xlim=c(0,15), ylim=c(0,15), lwd = 1,col="black", bg="orange"),ylab="JB1 Zfp57-/- log2(normalized read counts)", xlab="JB1 WT log2(normalized read counts)", bty = 'n') #brownish
lines(x = c(0,15), y = c(0,15), col = "#464647", lty = 2)
#Save manually as scatterund_bulk_mm10_t12_col_num_Zfp57dot_Inf2.svg
with(scatterund_re.clustimpgenes, text(WT,ZFP57KO, labels = Genes, pos = 3, cex=0.5, font=6))
#Save manually as idea for gene labels GENE_LABEL_scatterund_bulk_mm10_t12_col_num_Zfp57dot_Inf2.svg


scatterund_bulk_mm10_t12_column <- read.table("scatterund_bulk_mm10_t12_column.txt",header=F)
head(scatterund_bulk_mm10_t12_column)
dim(scatterund_bulk_mm10_t12_column)

#Unique_genes_removed duplicates due to Zfp57peaks overlap
sort -k1,1 -u scatterund_bulk_mm10_t12_column.txt > scatterund_bulk_mm10_t12_column_uniq.txt
scatterund_bulk_mm10_t12_column_uniq <- read.table("scatterund_bulk_mm10_t12_column_uniq.txt",header=F)
head(scatterund_bulk_mm10_t12_column_uniq)
dim(scatterund_bulk_mm10_t12_column_uniq)

fgrep -f scatterund_bulk_mm10_t12_col_num.Zfp57dot_clustimpgenes.txt scatterund_re.txt -w | awk '{print $2"\t"$3"\t"$4}' > scatterund_re.clustimpgenes.txt
scatterund_re.clustimpgenes <- read.table("scatterund_re.clustimpgenes.txt",header=F)
colnames(scatterund_re.clustimpgenes) <- c("Genes", "WT", "ZFP57KO")
head(scatterund_re.clustimpgenes)
dim(scatterund_re.clustimpgenes)
#Without Log#Not required Because data become skewed to nearby origin points
awk '{print $3"\t"$4"\t"$5"\t"$2"\t"$1"\t"$6"\t"$7}' ddsNormcountsavg_gene_WT_KO.pre.txt | grep ENS | sort -k1,1 -k2,2n > all_ddsNormcountsavg_gene_WT_KO.sorted.txt

#188 imprinted_gene_name.txt
fgrep -f imprinted_gene_name.txt all_ddsNormcountsavg_gene_WT_KO.sorted.txt -w > imprinted_ddsNormcountsavg_gene_WT_KO.sorted.txt 
fgrep -f imprinted_gene_name.txt all_ddsNormcountsavg_gene_WT_KO.sorted.txt -w -v > Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.txt

#upregulated
fgrep -f deseq2_results_res0.05_sorted.upregulated.ensid all_ddsNormcountsavg_gene_WT_KO.sorted.txt -w > all_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.txt
fgrep -f deseq2_results_res0.05_sorted.upregulated.ensid imprinted_ddsNormcountsavg_gene_WT_KO.sorted.txt -w > imprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.txt
fgrep -f deseq2_results_res0.05_sorted.upregulated.ensid Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.txt -w > Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.txt

#downregulated
fgrep -f deseq2_results_res0.05_sorted.downregulated.ensid all_ddsNormcountsavg_gene_WT_KO.sorted.txt -w > all_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.txt
fgrep -f deseq2_results_res0.05_sorted.downregulated.ensid imprinted_ddsNormcountsavg_gene_WT_KO.sorted.txt -w > imprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.txt
fgrep -f deseq2_results_res0.05_sorted.downregulated.ensid Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.txt -w > Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.txt


#Non-deregulated
fgrep -f deseq2_results_res0.05_sorted.ensid all_ddsNormcountsavg_gene_WT_KO.sorted.txt -w -v > all_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.txt
fgrep -f deseq2_results_res0.05_sorted.ensid imprinted_ddsNormcountsavg_gene_WT_KO.sorted.txt -w -v > imprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.txt
fgrep -f deseq2_results_res0.05_sorted.ensid Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.txt -w -v > Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.txt


##Downregulated
#All
bedtools closest -a all_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""downregulatedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > all_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a all_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""downregulatedNotunder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > all_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
#Imprinted: 
bedtools closest -a imprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""downregulatedImprintedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > imprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a imprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""downregulatedImprintedNotUnder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > imprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
#Nonimprinted: 
bedtools closest -a Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""downregulatedNotImprintedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""downregulatedNotImprintedNotUnder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt

##Upregulated
#All
bedtools closest -a all_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""upregulatedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > all_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a all_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""upregulatedNotunder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > all_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
#Imprinted: 
bedtools closest -a imprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""upregulatedImprintedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > imprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a imprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""upregulatedImprintedNotUnder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > imprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
#Nonimprinted: 
bedtools closest -a Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""upregulatedNotImprintedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""upregulatedNotImprintedNotUnder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt

#Nonderegulated
#All
bedtools closest -a all_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""NonDeregulatedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > all_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a all_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""NonDeregulatedNotunder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > all_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
#Imprinted: 
bedtools closest -a imprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""NonDeregulatedImprintedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > imprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a imprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""NonDeregulatedImprintedNotUnder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > imprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
#Nonimprinted: 
bedtools closest -a Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""NonDeregulatedNotImprintedUnder100Kb"}' | sort -k9,9n | awk '$9 <= 100000 {print $0}' > Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
bedtools closest -a Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$13"\t""NonDeregulatedNotImprintedNotUnder100Kb"}' | sort -k9,9n | awk '$9 > 100000 {print $0}' > Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt

wc -l all_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l all_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l imprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l imprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l all_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l all_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l imprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l imprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l all_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l all_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l imprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l imprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt
wc -l Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt

cat imprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt imprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.upregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt imprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt imprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.downregulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt  imprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt imprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt  Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt Nonimprinted_ddsNormcountsavg_gene_WT_KO.sorted.NOTderegulated.Notunder100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt | sort -k5,5 > scatterUNDbulk_mm10_t12.txt

awk '{print $4"\t"$6"\t"$7"\t"$10}'  scatterUNDbulk_mm10_t12.txt > scatterUND_bulk_mm10_t12_column.txt
awk  '{print $4}' scatterUND_bulk_mm10_t12_column.txt | sort -k1,1 -u #print all possibilities

downregulatedImprintedNotUnder100Kb 1
downregulatedImprintedUnder100Kb 2
downregulatedNotImprintedNotUnder100Kb 3
downregulatedNotImprintedUnder100Kb 4
NonDeregulatedImprintedNotUnder100Kb 9
NonDeregulatedImprintedUnder100Kb  10
NonDeregulatedNotImprintedNotUnder100Kb 11
NonDeregulatedNotImprintedUnder100Kb  12
upregulatedImprintedNotUnder100Kb  5
upregulatedImprintedUnder100Kb 6
upregulatedNotImprintedNotUnder100Kb  7
upregulatedNotImprintedUnder100Kb 8
Zfp57 13

#Replace effect name with numbers and remove -inf containing genes and save as  scatterUND_bulk_mm10_t12_col_num.txt #-inf not found
sed 's/downregulatedImprintedNotUnder100Kb/1/g' scatterUND_bulk_mm10_t12_column.txt | sed 's/downregulatedImprintedUnder100Kb/2/g' | sed 's/downregulatedNotImprintedNotUnder100Kb/3/g' | sed 's/downregulatedNotImprintedUnder100Kb/4/g' | sed 's/upregulatedImprintedNotUnder100Kb/5/g' | sed 's/upregulatedImprintedUnder100Kb/6/g' | sed 's/upregulatedNotImprintedNotUnder100Kb/7/g' | sed 's/upregulatedNotImprintedUnder100Kb/8/g' | sed 's/NonDeregulatedImprintedNotUnder100Kb/9/g' | sed 's/NonDeregulatedImprintedUnder100Kb/10/g' | sed 's/NonDeregulatedNotImprintedNotUnder100Kb/11/g' | sed 's/NonDeregulatedNotImprintedUnder100Kb/12/g' > scatterUND_bulk_mm10_t12_col_num.txt

#Assign ZFP57
grep Zfp57 scatterUND_bulk_mm10_t12_column.txt -w   > ZFP57DOT.txt
#Zfp57	336.033812598925	146.716091297545	downregulatedNotImprintedNotUnder100Kb
#vim ZFP57DOT.txt: downregulatedNotImprintedNotUnder100Kb to 13 and Zfp57 to ZFP57DOT
cat scatterUND_bulk_mm10_t12_col_num.txt ZFP57DOT.txt > scatterUND_bulk_mm10_t12_col_num.ZFP57DOT.txt
grep INF2 scatterUND_bulk_mm10_t12_col_num.ZFP57DOT.txt > INF2.txt
grep Inf scatterUND_bulk_mm10_t12_col_num.ZFP57DOT.txt -v > scatterUND_bulk_mm10_t12_col_num.ZFP57DOT_Inf.txt
cat scatterUND_bulk_mm10_t12_col_num.ZFP57DOT_Inf.txt INF2.txt > scatterUND_bulk_mm10_t12_col_num.ZFP57DOT_INF2.txt

library(ggplot2)
library(ggrepel)
scatterUND <- read.table("scatterUND_bulk_mm10_t12_col_num.ZFP57DOT_INF2.txt",header=F)
head(scatterUND)
scatterUND_re <- scatterUND
rownames(scatterUND_re)
colnames(scatterUND_re)
head(scatterUND_re)
dim(scatterUND_re)
write.table(scatterUND_re, "scatterUND_re.txt", quote = F, append = F)
scatterUND1 = scatterUND_re
head(scatterUND1)
dim(scatterUND1)
colnames(scatterUND1) <- c("Gene", "WT","ZFP57KO","Color")
head(scatterUND1)
tail(scatterUND1)
dim(scatterUND1)
scatterUND1 <- data.frame(scatterUND1)
with(scatterUND1, plot(WT, ZFP57KO, pch=21, cex = 0.2, main="Scatter Plot", xlim=c(0,10000), ylim=c(0,10000),lwd = 0.3,col="#525254", bg="grey",ylab="JB1 Zfp57-/- (normalized read counts)", xlab="JB1 WT (normalized read counts)", bty = 'n'))
with(subset(scatterUND1, Color==7), points(WT, ZFP57KO, pch=21, cex = 0.4, xlim=c(0,10000), ylim=c(0,10000), lwd = 0.5,col="#9fc5e8ff", bg="#9fc5e8ff"),ylab="", xlab="", bty = 'n') #reddish
with(subset(scatterUND1, Color==8), points(WT, ZFP57KO, pch=21, cex = 0.4, xlim=c(0,10000), ylim=c(0,10000), lwd = 0.5,col="#9fc5e8ff", bg="#9fc5e8ff"),ylab="", xlab="", bty = 'n') #yellowish
with(subset(scatterUND1, Color==5), points(WT, ZFP57KO, pch=21, cex = 0.6, xlim=c(0,10000), ylim=c(0,10000), lwd = 1,col="#0b5394ff", bg="#0b5394ff"),ylab="", xlab="", bty = 'n') #greenish
with(subset(scatterUND1, Color==6), points(WT, ZFP57KO, pch=21, cex = 0.6, xlim=c(0,10000), ylim=c(0,10000), lwd = 1,col="#0b5394ff", bg="#0b5394ff"),ylab="", xlab="", bty = 'n') #bluish
with(subset(scatterUND1, Color==3), points(WT, ZFP57KO, pch=21, cex = 0.4, xlim=c(0,10000), ylim=c(0,10000), lwd = 0.5,col="#f4ccccff", bg="#f4ccccff"),ylab="", xlab="", bty = 'n') #reddish
with(subset(scatterUND1, Color==4), points(WT, ZFP57KO, pch=21, cex = 0.4, xlim=c(0,10000), ylim=c(0,10000), lwd = 0.5,col="#f4ccccff", bg="#f4ccccff"),ylab="", xlab="", bty = 'n') #yellowish
with(subset(scatterUND1, Color==1), points(WT, ZFP57KO, pch=21, cex = 0.6, xlim=c(0,10000), ylim=c(0,10000), lwd = 1,col="#cc0000ff", bg="#cc0000ff"),ylab="", xlab="", bty = 'n') #greenish
with(subset(scatterUND1, Color==2), points(WT, ZFP57KO, pch=21, cex = 0.6, xlim=c(0,10000), ylim=c(0,10000), lwd = 1,col="#cc0000ff", bg="#cc0000ff"),ylab="", xlab="", bty = 'n') #bluish

with(subset(scatterUND1, Color==13), points(WT, ZFP57KO, pch=21, cex = 1.2, xlim=c(0,10000), ylim=c(0,10000), lwd = 1,col="black", bg="orange"),ylab="JB1 Zfp57-/- log2(normalized read counts)", xlab="JB1 WT log2(normalized read counts)", bty = 'n') #brownish
lines(x = c(0,10000), y = c(0,10000), col = "#464647", lty = 2)
#Save manually as scatterUND_bulk_mm10_t12_col_num.svg


# Functional analysis

#Note: IMPORT FILES ALREADY GENERATED AND RUN gostplot do not export new files
#A list of genes can be used to obtain a first biological insight by looking at some known annotation databases available online.
#By using the **gprofiler** package we can query some of these databases and looking at the cellular mechanisms mostly affected by the obtained gene list.
library(gProfileR)
library(gprofiler2)
gene_annotation <- read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20_rearranged.txt", header = FALSE)
head(gene_annotation)
colnames(gene_annotation) <- c("chr","start","end","id","gene")
gene_annotation <- data.frame(gene_annotation)
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2")
head(deseq2_results_res0.05_gene_rearranged)
#Upregulated and Downregulated
deseq2_results_res0.05_gene_rearranged <- read.table("deseq2_results_res0.05_gene_rearranged.txt", header = F)
head(deseq2_results_res0.05_gene_rearranged)
colnames(deseq2_results_res0.05_gene_rearranged) <- c("chr", "start", "end", "ensID", "Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "threshold")
head(deseq2_results_res0.05_gene_rearranged)
dim(deseq2_results_res0.05_gene_rearranged)
#Dr Claudia said no for this Filter Basemean>1st quartile summary(deseq2_results_res0.05_gene_rearranged$baseMean)
#deseq2_results_res0.05_gene_rearranged_baseMfilt <- deseq2_results_res0.05_gene_rearranged[which(deseq2_results_res0.05_gene_rearranged$baseMean > 54.31),]
#dim(deseq2_results_res0.05_gene_rearranged_baseMfilt)
#head(deseq2_results_res0.05_gene_rearranged_baseMfilt)
#write.table(deseq2_results_res0.05_gene_rearranged_baseMfilt$Gene,"deseq2_results_res0.05_gene_rearranged_baseMfilt.gene",quote = F, append = F, col.names = F, row.names = F)

#File preparation for gprofiler
#diff_0.5.up <- deseq2_results_res0.05_gene_rearranged_baseMfilt[which(deseq2_results_res0.05_gene_rearranged_baseMfilt$log2FoldChange > 0),]
diff_0.5.up <- deseq2_results_res0.05_gene_rearranged[which(deseq2_results_res0.05_gene_rearranged$log2FoldChange > 0),]
diff_0.5.up.id <- data.frame(diff_0.5.up$Gene)
colnames(diff_0.5.up.id) <- "id"
head(diff_0.5.up.id)
diff_0.5.up.symbol <- as.character(diff_0.5.up.id$id)
write.table(diff_0.5.up.symbol, "diff_0.5.up.symbol.txt", row.names = F, sep = "\t", quote = F, append = F)
diff_0.5.down <- deseq2_results_res0.05_gene_rearranged[which(deseq2_results_res0.05_gene_rearranged$log2FoldChange < 0),]
diff_0.5.down.id <- data.frame(diff_0.5.down$Gene)
colnames(diff_0.5.down.id) <- "id"
diff_0.5.down.symbol <- as.character(diff_0.5.down.id$id)
write.table(diff_0.5.down.symbol, "diff_0.5.down.symbol.txt", row.names = F, sep = "\t", quote = F, append = F)
library(gprofiler2)
gost.diff_0.5.up <- gost(query=diff_0.5.up.symbol, organism="mmusculus")
gostplot(gost.diff_0.5.up)
gost.diff_0.5.down <- gost(query=diff_0.5.down.symbol, organism="mmusculus")
gostplot(gost.diff_0.5.down)

#Online version also ran and results were comparable see csv files

############----------------------
#Note: p_value output of gost function is Hypergeometric p-value after correction for multiple testing, see https://biit.cs.ut.ee/gprofiler_beta/page/apis#gost_query_results
#Sort by p-value, Use PRE-RAN FILES
#Up
#gost.diff_0.5.up.res.sorted <- read.delim("gost.diff_0.5.up.res.sorted.txt")
#gost.diff_0.5.down.res.sorted <- read.delim("gost.diff_0.5.down.res.sorted.txt")
gost.diff_0.5.up.res.sorted <- gost.diff_0.5.up$result[order(gost.diff_0.5.up$result$p_value),]
head(gost.diff_0.5.up.res.sorted)
write.table(as.matrix(gost.diff_0.5.up.res.sorted),"gost.diff_0.5.up.res.sorted.txt",quote = F, append = F, sep = "\t")
write.xlsx(as.matrix(gost.diff_0.5.up.res.sorted),"gost.diff_0.5.up.res.sorted.xlsx")

gost.diff_0.5.down.res.sorted <- gost.diff_0.5.down$result[order(gost.diff_0.5.down$result$p_value),]
head(gost.diff_0.5.down.res.sorted)
write.table(as.matrix(gost.diff_0.5.down.res.sorted),"gost.diff_0.5.down.res.sorted.txt",quote = F, append = F, sep = "\t")
#GO:BP
gost.diff_0.5.up.res.sorted_gobp <- gost.diff_0.5.up.res.sorted[which(gost.diff_0.5.up.res.sorted$source == "GO:BP"),]
head(gost.diff_0.5.up.res.sorted_gobp)
dim(gost.diff_0.5.up.res.sorted_gobp)
gost.diff_0.5.up.res.sorted_gobp_bar <- gost.diff_0.5.up.res.sorted_gobp[,c(11,3)]
head(gost.diff_0.5.up.res.sorted_gobp_bar)
gost.diff_0.5.up.res.sorted_gobp_bar_top <- head(gost.diff_0.5.up.res.sorted_gobp_bar,10)
gost.diff_0.5.up.res.sorted_gobp_bar_top$term_name <- gsub(' ', '.', gost.diff_0.5.up.res.sorted_gobp_bar_top$term_name)
gost.diff_0.5.up.res.sorted_gobp_bar_top$p_value <- -log10(gost.diff_0.5.up.res.sorted_gobp_bar_top$p_value)
gost.diff_0.5.up.res.sorted_gobp_bar_top <- data.frame(gost.diff_0.5.up.res.sorted_gobp_bar_top)

library(ggpubr)
ggbarplot(gost.diff_0.5.up.res.sorted_gobp_bar_top, 
          x = "term_name", 
          y = "p_value", 
          color = "#5d93c4ff", 
          fill = "#5d93c4ff" ,
          sort.by.groups = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "-log10(p.value)",
          xlab = "Pathways",
          legend.title = "Pathways",
          lab.size = 9,
          sort.val = "asc",
          rotate = TRUE, 
          position = position_dodge(),
          ggtheme = theme_bw())

ggsave("gost.diff_0.5.up.res.sorted_gobp_bar_top.svg", width=12, height=8, units="cm", dpi=96)

#Down

#GO:BP
gost.diff_0.5.down.res.sorted_gobp <- gost.diff_0.5.down.res.sorted[which(gost.diff_0.5.down.res.sorted$source == "GO:BP"),]
head(gost.diff_0.5.down.res.sorted_gobp)
dim(gost.diff_0.5.down.res.sorted_gobp)
gost.diff_0.5.down.res.sorted_gobp_bar <- gost.diff_0.5.down.res.sorted_gobp[,c(11,3)]
head(gost.diff_0.5.down.res.sorted_gobp_bar)
gost.diff_0.5.down.res.sorted_gobp_bar_top <- head(gost.diff_0.5.down.res.sorted_gobp_bar,10)
gost.diff_0.5.down.res.sorted_gobp_bar_top$term_name <- gsub(' ', '.', gost.diff_0.5.down.res.sorted_gobp_bar_top$term_name)
gost.diff_0.5.down.res.sorted_gobp_bar_top$p_value <- -log10(gost.diff_0.5.down.res.sorted_gobp_bar_top$p_value)
gost.diff_0.5.down.res.sorted_gobp_bar_top <- data.frame(gost.diff_0.5.down.res.sorted_gobp_bar_top)

library(ggpubr)
ggbarplot(gost.diff_0.5.down.res.sorted_gobp_bar_top, 
          x = "term_name", 
          y = "p_value", 
          color = "#e58e8eff", 
          fill = "#e58e8eff" ,
          sort.by.grodowns = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "-log10(p.value)",
          xlab = "Pathways",
          legend.title = "Pathways",
          lab.size = 9,
          rotate = TRUE, position = position_dodge(),
          sort.val = "asc",
          ggtheme = theme_bw())

ggsave("gost.diff_0.5.down.res.sorted_gobp_bar_top.svg", width=13.5, height=8, units="cm", dpi=96)

#Downregulated under 100Kb
#all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt (188 features) some genes are present more times due to Zfp57 peaks overlap
#Filterout mitochondrial chromosome genes and Y (177 leftover)
grep chrM all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10.txt -v > all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10_subchrM.txt
grep chrY all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10_subchrM.txt -v > all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10_subchrMY.txt

fgrep -f deseq2_results_res0.05_gene_rearranged.gene all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10_subchrMY.txt -w > all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10_subchrMY_filt.txt
#gProfiler2
diff_0.5.downun100kb <- read.table("all_ddsNormcountsavgLog_gene_WT_KO.sorted.downregulated.under100Kb.JB1_Zfp57_overlapped_Kap1_mm10_subchrMY_filt.txt", header =F)
head(diff_0.5.downun100kb)
colnames(diff_0.5.downun100kb) <- c("chr","start","end","Gene","id","WT","ZFP57KO","Dot","Distance","Effect")
diff_0.5.downun100kb.id <- data.frame(diff_0.5.downun100kb$Gene)
colnames(diff_0.5.downun100kb.id) <- "id"
diff_0.5.downun100kb.symbol <- unique(as.character(diff_0.5.downun100kb.id$id))#deduplicate symbol
write.table(diff_0.5.downun100kb.symbol, "diff_0.5.downun100kb.symbol.txt", quote=F, append=F)
#Down under 100Kb
gost.diff_0.5.downun100kb <- gost(query=diff_0.5.downun100kb.symbol, organism="mmusculus")
gostplot(gost.diff_0.5.downun100kb)
#Sort by p-value, Use pre-ran files
#gost.diff_0.5.downun100kb.res.sorted <- read.delim("gost.diff_0.5.downun100kb.res.sorted.txt")
gost.diff_0.5.downun100kb.res.sorted <- gost.diff_0.5.downun100kb$result[order(gost.diff_0.5.downun100kb$result$p_value),]
head(gost.diff_0.5.downun100kb.res.sorted)
#write.table(as.matrix(gost.diff_0.5.downun100kb.res.sorted),"gost.diff_0.5.downun100kb.res.sorted.txt",quote = F, append = F, sep = "\t")


#GO:BP
gost.diff_0.5.downun100kb.res.sorted_gobp <- gost.diff_0.5.downun100kb.res.sorted[which(gost.diff_0.5.downun100kb.res.sorted$source == "GO:BP"),]
head(gost.diff_0.5.downun100kb.res.sorted_gobp)
dim(gost.diff_0.5.downun100kb.res.sorted_gobp)
gost.diff_0.5.downun100kb.res.sorted_gobp_bar <- gost.diff_0.5.downun100kb.res.sorted_gobp[,c(11,3)]
head(gost.diff_0.5.downun100kb.res.sorted_gobp_bar)
gost.diff_0.5.downun100kb.res.sorted_gobp_bar_top <- head(gost.diff_0.5.downun100kb.res.sorted_gobp_bar,20)
gost.diff_0.5.downun100kb.res.sorted_gobp_bar_top$term_name <- gsub(' ', '.', gost.diff_0.5.downun100kb.res.sorted_gobp_bar_top$term_name)
gost.diff_0.5.downun100kb.res.sorted_gobp_bar_top$p_value <- -log10(gost.diff_0.5.downun100kb.res.sorted_gobp_bar_top$p_value)
gost.diff_0.5.downun100kb.res.sorted_gobp_bar_top <- data.frame(gost.diff_0.5.downun100kb.res.sorted_gobp_bar_top)

library(ggpubr)
ggbarplot(gost.diff_0.5.downun100kb.res.sorted_gobp_bar_top, 
          x = "term_name", 
          y = "p_value", 
          color = "#e58e8eff", 
          fill = "#e58e8eff" ,
          sort.by.grodowns = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "-log10(p.value)",
          xlab = "Pathways",
          legend.title = "Pathways",
          lab.size = 9,
          rotate = TRUE, position = position_dodge(),
          sort.val = "asc",
          ggtheme = theme_bw())

ggsave("gost.diff_0.5.downun100kb.res.sorted_gobp_bar_top.svg", width=13.5, height=8, units="cm", dpi=96)

#extract details of down regulated gene under 100Kb (for excel files)
head(diff_0.5.downun100kb)
head(deseq2_results_res0.05_gene_rearranged)
diff_0.5.downun100kb.deg <- merge(diff_0.5.downun100kb,deseq2_results_res0.05_gene_rearranged, by ="Gene", all.x=FALSE,sort = F)
dim(diff_0.5.downun100kb.deg)
head(diff_0.5.downun100kb.deg)
diff_0.5.downun100kb.deg.re <- diff_0.5.downun100kb.deg[,c(11:20, 1:10)]
#There are some genes which were duplicated while intersecting with ZFP57/KAP1 peaks, keep it as such and delet only in excel sheet
library(WriteXLS) or use library("writexl")write_xlsx()

write_xlsx(diff_0.5.downun100kb.deg.re, "diff_0.5.downun100kb.deg.xlsx")
#Take deduplicated set of these genes 138 from excel (use venny) and paste in gprofiler onlinde version with mmusculus
#Check if Both online and R version give same results #If same go ahead
#Click Detailed Results
#Click CSV
#Copy downloaded csv file in current folder
diff_0.5.downun100kb.deg.re.gPro <- read.csv("gProfiler_mmusculus_1-23-2021_4-09-02 PM__intersections.csv")
head(diff_0.5.downun100kb.deg.re.gPro)
diff_0.5.downun100kb.deg.re.gPro_gobp <- diff_0.5.downun100kb.deg.re.gPro[which(diff_0.5.downun100kb.deg.re.gPro$source == "GO:BP"),]
head(diff_0.5.downun100kb.deg.re.gPro_gobp)
dim(diff_0.5.downun100kb.deg.re.gPro_gobp)
#Take 20 candidates and term_id and intersections
diff_0.5.downun100kb.deg.re.gPro_gobp_re <- diff_0.5.downun100kb.deg.re.gPro_gobp[c(1:20),c(3,10)]
head(diff_0.5.downun100kb.deg.re.gPro_gobp_re)
dim(diff_0.5.downun100kb.deg.re.gPro_gobp_re)
write_xlsx(diff_0.5.downun100kb.deg.re.gPro_gobp_re,"diff_0.5.downun100kb.deg.re.gPro_gobp_re.xlsx")

library(splitstackshape)
diff_0.5.downun100kb.deg.re.gPro_gobp_resep <- data.frame(cSplit(diff_0.5.downun100kb.deg.re.gPro_gobp_re, "intersections", ","))
head(diff_0.5.downun100kb.deg.re.gPro_gobp_resep,1)
dim(diff_0.5.downun100kb.deg.re.gPro_gobp_resep)
rownames(diff_0.5.downun100kb.deg.re.gPro_gobp_resep) <- diff_0.5.downun100kb.deg.re.gPro_gobp_resep[,1]
head(diff_0.5.downun100kb.deg.re.gPro_gobp_resep,1)
#Checked manually till here alright!
library(reshape2)
diff_0.5.downun100kb.deg.re.gPro_gobp_resep1 <- melt(diff_0.5.downun100kb.deg.re.gPro_gobp_resep, id.vars="term_id") #When you melt, you are combining multiple columns into one.
#Warning message:attributes are not identical across measure variables; they will be dropped 
head(diff_0.5.downun100kb.deg.re.gPro_gobp_resep1)
dim(diff_0.5.downun100kb.deg.re.gPro_gobp_resep1) #20 *72 =1440, so melt works in separating and merging
#Checked manully using excel
#Rearrange column
diff_0.5.downun100kb.deg.re.gPro_gobp_resep2 <- diff_0.5.downun100kb.deg.re.gPro_gobp_resep1[,c(3,1,2)]
head(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2)
colnames(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2) <- c("Genes", "term_id","intersections")
head(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2)
#Sort by Gene alphabets
diff_0.5.downun100kb.deg.re.gPro_gobp_resep2sort <- diff_0.5.downun100kb.deg.re.gPro_gobp_resep2[order(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2$Genes),]
head(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2sort)
diff_0.5.downun100kb.deg.re.gPro_gobp_resep2sort <- data.frame(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2sort)
write.table(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2sort, "diff_0.5.downun100kb.deg.re.gPro_gobp_resep2sort.txt", sep = "\t", quote = F, append = F, row.names = F)
#Read throught the output
diff_0.5.downun100kb.deg.re.gPro_gobp_resep2sort <- read.table("diff_0.5.downun100kb.deg.re.gPro_gobp_resep2sort.txt", header = T)

diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab <- table(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2sort$Genes, diff_0.5.downun100kb.deg.re.gPro_gobp_resep2sort$term_id)
head(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab)

diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1 <- as.data.frame.matrix(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab)
head(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1)
tdiff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1 <- t(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1)
tdiff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1 <- data.frame(tdiff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1)
head(tdiff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1)
tdiff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1["term_id"] <- rownames(tdiff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1)
head(tdiff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1)
dim(tdiff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1)#There are 98 genes and 20 terms in the list Downreg100kb_functional_analysis_manaully checked.xlsx
#dim gives 20 99 as one is term_id

#Go back to you GO:BP term list
head(diff_0.5.downun100kb.deg.re.gPro_gobp)
dim(diff_0.5.downun100kb.deg.re.gPro_gobp)#They are 40
#Take top 20 GO:BP (row 1:20 data already sorted by adj p value) and take , termname, termid,-log10(padj value), intersections 2,3,5,10
diff_0.5.downun100kb.deg.re.gPro_gobp_re2 <- diff_0.5.downun100kb.deg.re.gPro_gobp[c(1:20),c(2:3,5,10)]
head(diff_0.5.downun100kb.deg.re.gPro_gobp_re2)
write_xlsx(diff_0.5.downun100kb.deg.re.gPro_gobp_re2,"diff_0.5.downun100kb.deg.re.gPro_gobp_re2.xlsx")

#Combined my binary counts and termname by using term_id
tdiff_0.5.downun100kb.deg.re.gPro_gobp_re2.deg <- merge(diff_0.5.downun100kb.deg.re.gPro_gobp_re2,tdiff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1, by ="term_id", all.x=T,sort = F)
dim(tdiff_0.5.downun100kb.deg.re.gPro_gobp_re2.deg)
head(tdiff_0.5.downun100kb.deg.re.gPro_gobp_re2.deg)
write_xlsx(tdiff_0.5.downun100kb.deg.re.gPro_gobp_re2.deg,"tdiff_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.xlsx")
dim(tdiff_0.5.downun100kb.deg.re.gPro_gobp_re2.deg)
#Take term_id and all 98 genes binary counts
tdiff_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.re <- tdiff_0.5.downun100kb.deg.re.gPro_gobp_re2.deg[,c(1,5:length(tdiff_0.5.downun100kb.deg.re.gPro_gobp_re2.deg))]
dim(tdiff_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.re)
head(tdiff_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.re)
rownames(tdiff_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.re) <- tdiff_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.re$term_id
tdiff_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.re<- tdiff_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.re[,-1]
#Transpose back the gene and term_id
diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2 <- data.frame(t(tdiff_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.re))
head(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2)
diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2["uppercaseGene"] <- rownames(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2)
head(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2)
dim(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2)#98 21, 98 genes, 20 GO:BP terms and 1 uppercaseGene column
#Rearrange
diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2 <- diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2[,c(21,1:20)]
diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2 <- data.frame(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2)
write_xlsx(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2,"diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2.xlsx")

#Now Take deregulated genes undwer 100 kb from Zfp57/KAP1 peaks dataset
head(diff_0.5.downun100kb.deg.re)
dim(diff_0.5.downun100kb.deg.re)
diff_0.5.downun100kb.deg.re.withuppercase <- diff_0.5.downun100kb.deg.re
diff_0.5.downun100kb.deg.re.withuppercase["uppercaseGene"] <- data.frame((toupper(diff_0.5.downun100kb.deg.re$Gene)))
head(diff_0.5.downun100kb.deg.re.withuppercase)
write.table(diff_0.5.downun100kb.deg.re.withuppercase$uppercaseGene, "diff_0.5.downun100kb.deg.re.withuppercase.id",quote = F,append = F, row.names = F, col.names = F )
#Merge the 98 genes, 20 GO:BP terms and 1 uppercaseGene column with Zfp57/Kap1 100kb information
diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1.deg <- merge(diff_0.5.downun100kb.deg.re.withuppercase,diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2, by ="uppercaseGene", all.x=T,sort = F)
dim(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1.deg)
head(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1.deg)
#The informartion for all genes will not be available so NA will be depicted
#manualy checked see Downreg100kb_functional_analysis_manaullychecked.ods
write_xlsx(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1.deg,"diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1.deg.xlsx")
write.table(diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1.deg,"diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1.deg.txt",row.names = F, quote=F, append = F, sep = "\t")

#Label G:BP terms related to required features

#Take all candidates and term_id and intersections
diff_tot_0.5.downun100kb.deg.re.gPro_gobp_re <- diff_0.5.downun100kb.deg.re.gPro_gobp[,c(3,10)]
head(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_re)
dim(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_re)
write_xlsx(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_re,"diff_tot_0.5.downun100kb.deg.re.gPro_gobp_re.xlsx")

library(splitstackshape)
diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep <- data.frame(cSplit(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_re, "intersections", ","))
head(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep,1)
dim(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep)
rownames(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep) <- diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep[,1]
head(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep,1)
#Checked manually till here alright!
library(reshape2)
diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep1 <- melt(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep, id.vars="term_id") #When you melt, you are combining multiple columns into one.
#Warning message:attributes are not identical across measure variables; they will be dropped 
head(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep1)
dim(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep1) #40 *72 =2880, so melt works in separating and merging
#Checked manully using excel
#Rearrange column
diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2 <- diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep1[,c(3,1,2)]
head(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2)
colnames(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2) <- c("Genes", "term_id","intersections")
head(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2)
#Sort by Gene alphabets
diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2sort <- diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2[order(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2$Genes),]
head(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2sort)
diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2sort <- data.frame(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2sort)
write.table(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2sort, "diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2sort.txt", sep = "\t", quote = F, append = F, row.names = F)
#Read throught the output
diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2sort <- read.table("diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2sort.txt", header = T)

diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab <- table(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2sort$Genes, diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2sort$term_id)
head(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab)

diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1 <- as.data.frame.matrix(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab)
head(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1)
tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1 <- t(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1)
tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1 <- data.frame(tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1)
head(tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1)
tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1["term_id"] <- rownames(tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1)
head(tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1)
dim(tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1)#There are 98 genes and 40 terms in the list Downreg100kb_functional_analysis_manaully checked.xlsx
#dim gives 20 99 as one is term_id

#Go back to you GO:BP term list
head(diff_0.5.downun100kb.deg.re.gPro_gobp)
dim(diff_0.5.downun100kb.deg.re.gPro_gobp)#They are 40
#Take all GO:BP (data already sorted by adj p value) and take , termname, termid,-log10(padj value), intersections 2,3,5,10
diff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2 <- diff_0.5.downun100kb.deg.re.gPro_gobp[,c(2:3,5,10)]
head(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2)
dim(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2)
write_xlsx(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2,"diff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2.xlsx")

#Combined my binary counts and termname by using term_id
tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2.deg <- merge(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2,tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1, by ="term_id", all.x=T,sort = F)
dim(tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2.deg)
head(tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2.deg)
write_xlsx(tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2.deg,"tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.xlsx")
dim(tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2.deg)
#Take term_id and all 98 genes binary counts
tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.re <- tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2.deg[,c(1,5:length(tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2.deg))]
dim(tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.re)
head(tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.re)
rownames(tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.re) <- tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.re$term_id
tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.re<- tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.re[,-1]
#Transpose back the gene and term_id
diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2 <- data.frame(t(tdiff_tot_0.5.downun100kb.deg.re.gPro_gobp_re2.deg.re))
head(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2)
diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2["uppercaseGene"] <- rownames(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2)
head(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2)
dim(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2)#98 41, 98 genes, 40 GO:BP terms and 1 uppercaseGene column
#Rearrange
diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2 <- diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2[,c(41,1:40)]
diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2 <- data.frame(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2)
write_xlsx(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2,"diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2.xlsx")

#Now Take deregulated genes undwer 100 kb from Zfp57/KAP1 peaks dataset
head(diff_0.5.downun100kb.deg.re)
dim(diff_0.5.downun100kb.deg.re)
diff_tot_0.5.downun100kb.deg.re.withuppercase <- diff_0.5.downun100kb.deg.re
diff_tot_0.5.downun100kb.deg.re.withuppercase["uppercaseGene"] <- data.frame((toupper(diff_0.5.downun100kb.deg.re$Gene)))
head(diff_tot_0.5.downun100kb.deg.re.withuppercase)
write.table(diff_tot_0.5.downun100kb.deg.re.withuppercase$uppercaseGene, "diff_tot_0.5.downun100kb.deg.re.withuppercase.id",quote = F,append = F, row.names = F, col.names = F )
#Merge the 98 genes, 40 GO:BP terms and 1 uppercaseGene column with Zfp57/Kap1 100kb information
diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1.deg <- merge(diff_tot_0.5.downun100kb.deg.re.withuppercase,diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab2, by ="uppercaseGene", all.x=T,sort = F)
dim(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1.deg)
head(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1.deg)
#The informartion for all genes will not be available so NA will be depicted
#manualy checked see Downreg100kb_functional_analysis_manaullychecked.ods
write_xlsx(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1.deg,"diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1.deg.xlsx")
write.table(diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1.deg,"diff_tot_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1.deg.txt",row.names = F, quote=F, append = F, sep = "\t")

#Label GO:BP terms related to required features


#Open gprofiler output csv

#Get imprinted genes table
fgrep -f imprinted_gene_name.txt diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1.deg.txt -w > imprinted_diff_0.5.downun100kb.deg.re.gPro_gobp_resep2tab1.deg.txt

#Distance from closest ZFP57 peaks

#XXXXXXXXXXXXXX#####Do not run theis command for now, Do not filter Filter non-windows/biased regions : For detailed description see Oct2019 in documents folder 
#XXXXXXXXXXXXXXXXX######sort -k1,1 -k2,2n all_deregulated_genes_deseq2_0.05.txt | grep ENS > all_deregulated_genes_deseq2_0.05.sorted1.txt
#XXXXXXXXXXXXXXXXXXXX######bedtools intersect -wa -wb -a all_deregulated_genes_deseq2_0.05.sorted1.txt -b combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort_coverage_toremove.txt -v > all_deregulated_genes_deseq2_0.05.sorted.txt
#Currently run this,
#188 imprinted_gene_name.txt
#For interesection with bedtools: 
grep chr deseq2_results_res0.05_gene_rearranged.txt > all_deregulated_genes_deseq2_0.05.txt
sort -k1,1 -k2,2n all_deregulated_genes_deseq2_0.05.txt | grep ENS > all_deregulated_genes_deseq2_0.05.sorted.txt
fgrep -f imprinted_gene_name.txt all_deregulated_genes_deseq2_0.05.sorted.txt -w > imprinted_deregulated_genes_deseq2_0.05.sorted.txt 
fgrep -f imprinted_gene_name.txt all_deregulated_genes_deseq2_0.05.sorted.txt -w -v > Nonimprinted_deregulated_genes_deseq2_0.05.sorted.txt

#Overlap JB1_Zfp57 (544) and JB1_Kap1 (6986) peaks and extract only JB1_Zfp57 peaks (457)
#bedtools intersect -wa -wb -a Zfp57_GSM1931786_mm10_peaks.bed -b Kap1_GSM1931787_mm10_peaks.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > Zfp57_overlapped_Kap1_mm10_peaks.bed

###### Do not run this command. do not filter deletion/biased regions for now, bedtools intersect -wa -wb -a deseq2_results_res0.05_gene_rearranged.txt -b combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort_coverage_toremove.txt -v > deseq2_results_res0.05_gene_rearranged_filtered.txt (3989)
#Important do not remove duplicated genes because distance of genes may vary from one peak to another, sort low to high use sort -k9,9n
#All: 
bedtools closest -a all_deregulated_genes_deseq2_0.05.sorted.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$18}' | sort -k9,9n > all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt
#Some peaks of Zfp57 are present in split for the same gene example Plagl1, so number will rise from 4912 to 4921. But I  removed the Duplicates during counting as in ./distanceplot_awk.sh
#Imprinted: 
bedtools closest -a imprinted_deregulated_genes_deseq2_0.05.sorted.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$18}' | sort -k9,9n > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt
#Nonimprinted: 
bedtools closest -a Nonimprinted_deregulated_genes_deseq2_0.05.sorted.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11"\t"$18}' | sort -k9,9n > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt

#Remove chrM and chrY as the distance were -1
grep chrM all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -v > all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrM.txt
grep chrY all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrM.txt -v > all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrMY.txt
grep chrM imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -v > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrM.txt
grep chrY imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrM.txt -v > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrMY.txt
grep chrM Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -v > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrM.txt
grep chrY Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrM.txt -v > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrMY.txt
grep chrM res_sorted_gene_rearranged.txt -v > res_sorted_gene_rearranged_M.txt
grep chrY res_sorted_gene_rearranged_M.txt -v > res_sorted_gene_rearranged_MY.txt
grep chrM deseq2_results_res0.05_gene_rearranged.txt -v > deseq2_results_res0.05_gene_rearranged_M.txt
grep chrY deseq2_results_res0.05_gene_rearranged_M.txt -v > deseq2_results_res0.05_gene_rearranged_MY.txt

#Run this shell script 
./distanceplot_awkre.sh > distanceplot_valuesre.txt
# Line plot with multiple groups
#distancepeaks <- c('0','10','100','1000','10000','100000','1000000','10000000','100000000','1000000000')
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/")
distanceplotvalue <- read.table("distanceplot_valuesre.txt", header=FALSE)
head(distanceplotvalue)
plotdistance <- data.frame(c('0','10','100','1000','10000','100000','1000000','10000000','100000000','1000000000'))
#Denominator original (deduplicated) due to peaks overlap are as follows All: 4908(4899), imprinted 55(50), Non imprinted 4853(4849)
plotdistance["All"] <- as.numeric(as.character(distanceplotvalue[2:11,]))*100/4899
plotdistance["Imprinted"] <- as.numeric(as.character(distanceplotvalue[13:22,]))*100/50
plotdistance["Nonmprinted"] <- as.numeric(as.character(distanceplotvalue[24:33,]))*100/4849
colnames(plotdistance) <- c("Group", "All", "Imprinted", "Nonimprinted")
plotdistance
rownames(plotdistance) <- plotdistance$Group
plotdistance["Group"] <- NULL
plotdistance1 <- stack(plotdistance)
head(plotdistance1)
plotdistance1["Distance"] <- c(0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9)
colnames(plotdistance1) <- c("Genes", "Group", "Distance")
plotdistance1 <- plotdistance1[,c(2,3,1)]
head(plotdistance1)
write.table(plotdistance1,"/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/plotdistance1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)
d <- ggplot(data=plotdistance1, aes(x=Distance, y=Genes, group=Group)) +
  geom_line(aes(color=Group), size=1)+theme_classic()
d+scale_color_manual(values=c("black", "#00b050ff", "orange"))+ 
  scale_y_continuous(breaks = seq(0, 100, by = 10)) + 
  scale_x_continuous(breaks = seq(0, 100, by = 1))+
  labs(title="Distance plot", y ="% of deregulated genes in JB1", x = "Distance from JB1 ZFP57 peaks (bp)")+
  theme(axis.text.x = element_text(color="black", size=10, angle=315))
ggsave("distanceplot.svg", width=12, height=8, units="cm", dpi=96)

#As discussed by Basilia remove "All" Black curve
plotdistanceNoAll <- plotdistance1[-(1:10),]
disImpNonImp <- ggplot(data=plotdistanceNoAll, aes(x=Distance, y=Genes, group=Group)) +
  geom_line(aes(color=Group), size=1)+theme_classic()
disImpNonImp +scale_color_manual(values=c("#00b050ff", "orange"))+ 
  scale_y_continuous(breaks = seq(0, 100, by = 10)) + 
  scale_x_continuous(breaks = seq(0, 100, by = 1))+
  labs(title="Distance plot", y ="% of deregulated genes in JB1", x = "Distance from JB1 ZFP57 peaks (bp)")+
  theme(axis.text.x = element_text(color="black", size=10, angle=315))
ggsave("distanceplot_onlyImpNonImp.svg", width=12, height=8, units="cm", dpi=96)



#Cluster Dotplot of 40 imprinted genes. 
#Basilia provided the an excel sheet list_for_dotplot_allelic_ratio_e_ratio_bulk.xlsx
#I took gene list names (saved as basilia_genes_for_clusterdotplot_as_bulk.txt) to check counts
basilia_genes_for_clusterdotplot_as_bulk <- read.table("basilia_genes_for_clusterdotplot_as_bulk.txt", header = F)
colnames(basilia_genes_for_clusterdotplot_as_bulk) <- c("Genes")
head(ddsNormcounts_gene)
dim(ddsNormcounts_gene)
ddsNormcounts_chr_gene1_clusterdotplotgenes <- merge(basilia_genes_for_clusterdotplot_as_bulk,ddsNormcounts_gene, by ="Genes", all.x=FALSE,sort = F)
head(ddsNormcounts_chr_gene1_clusterdotplotgenes)
dim(ddsNormcounts_chr_gene1_clusterdotplotgenes)
library("writexl") #Remember it does not print rowname, and also I dont need rownames here in this particular  case
write_xlsx(ddsNormcounts_chr_gene1_clusterdotplotgenes,"ddsNormcounts_chr_gene1_clusterdotplotgenes.xlsx")

#First calculate average because as Basilia discussed rep1 can be compared to both rep1 and rep2 not to only one
ddsNormcounts_chr_gene1_clusterdotplotgenes["WT.avg"] <- (ddsNormcounts_chr_gene1_clusterdotplotgenes$JB1_WT_Rep1+ddsNormcounts_chr_gene1_clusterdotplotgenes$JB1_WT_Rep2)/2
ddsNormcounts_chr_gene1_clusterdotplotgenes["ZFP57KO.avg"] <- (ddsNormcounts_chr_gene1_clusterdotplotgenes$JB1_ZFP57_KO_Rep1+ddsNormcounts_chr_gene1_clusterdotplotgenes$JB1_ZFP57_KO_Rep2)/2
head(ddsNormcounts_chr_gene1_clusterdotplotgenes)

#calculate ko/(wt+ko)
ddsNormcounts_chr_gene1_clusterdotplotgenes["Ratio.avg"] <- ddsNormcounts_chr_gene1_clusterdotplotgenes$ZFP57KO.avg/(ddsNormcounts_chr_gene1_clusterdotplotgenes$WT.avg+ddsNormcounts_chr_gene1_clusterdotplotgenes$ZFP57KO.avg)
head(ddsNormcounts_chr_gene1_clusterdotplotgenes)
dim(ddsNormcounts_chr_gene1_clusterdotplotgenes)
#Separate replicates
cluster_imp_gene_indiv_bulk_norm_are <- ddsNormcounts_chr_gene1_clusterdotplotgenes[,c(1,12)]
head(cluster_imp_gene_indiv_bulk_norm_are)
dim(cluster_imp_gene_indiv_bulk_norm_are)
rownames(cluster_imp_gene_indiv_bulk_norm_are) <- cluster_imp_gene_indiv_bulk_norm_are[,1]

cluster_imp_gene_indiv_bulk_norm_sortre <- data.frame(cluster_imp_gene_indiv_bulk_norm_are$Ratio.avg)
head(cluster_imp_gene_indiv_bulk_norm_sortre)
rownames(cluster_imp_gene_indiv_bulk_norm_sortre) <- rownames(cluster_imp_gene_indiv_bulk_norm_are) 
cluster_imp_gene_indiv_bulk_norm_sortre
colnames(cluster_imp_gene_indiv_bulk_norm_sortre) <- "Ratio.avg"

#cbind.data.frame(cluster_imp_gene_indiv_bulk_norm_are,cluster_imp_gene_indiv_bulk_norm_sortre )
#checked the parallel while assigning rename
cluster_imp_gene_indiv_bulk_norm_sortrebarplot <- as.matrix(t(cluster_imp_gene_indiv_bulk_norm_sortre))
head(cluster_imp_gene_indiv_bulk_norm_sortrebarplot)
dim(cluster_imp_gene_indiv_bulk_norm_sortrebarplot)
cluster_imp_gene_indiv_bulk_norm_sortrebarplot1 <- stack(cluster_imp_gene_indiv_bulk_norm_sortrebarplot)
head(cluster_imp_gene_indiv_bulk_norm_sortrebarplot1)
colnames(cluster_imp_gene_indiv_bulk_norm_sortrebarplot1) <- c("Group", "Gene", "Ratio")
cluster_imp_gene_indiv_bulk_norm_sortrebarplot1 <- data.frame(cluster_imp_gene_indiv_bulk_norm_sortrebarplot1)
head(cluster_imp_gene_indiv_bulk_norm_sortrebarplot1)
dim(cluster_imp_gene_indiv_bulk_norm_sortrebarplot1)
cluster_imp_gene_indiv_bulk_norm_sortrebarplot1["Group.1"] <- data.frame(rep("KOWT_ratio", 40))
head(cluster_imp_gene_indiv_bulk_norm_sortrebarplot1)
dim(cluster_imp_gene_indiv_bulk_norm_sortrebarplot1)
str(cluster_imp_gene_indiv_bulk_norm_sortrebarplot1)
cluster_imp_gene_indiv_bulk_norm_sortrebarplot1_half <- cluster_imp_gene_indiv_bulk_norm_sortrebarplot1
cluster_imp_gene_indiv_bulk_norm_sortrebarplot1_half["Ratio"] <- as.numeric(cluster_imp_gene_indiv_bulk_norm_sortrebarplot1_half$Ratio) - 0.5
head(cluster_imp_gene_indiv_bulk_norm_sortrebarplot1_half)
library(ggplot2)
# Dotplot for Cluster
# Change the position
pcluster <-ggplot(cluster_imp_gene_indiv_bulk_norm_sortrebarplot1_half, aes(x=Gene, y=Ratio, fill= Group.1)) + geom_hline(yintercept = c(-0.17,0,0.17), colour = "grey", linetype=c("dashed","solid","dashed"))+
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.8))
pcluster + scale_fill_manual(values=c("blue"))+ylim(c(-0.5,0.5))+
  scale_color_manual(values=c("black","black","black","black")) + theme_classic()

ggsave("dotplot_cluster_imp_gene_indiv_bulk_norm_sortrebarplot1_half.svg", width=100, height=10, units="cm", dpi=96)
ggsave("dotplot_cluster_imp_gene_indiv_bulk_norm_sortrebarplot1_half.jpg", width=100, height=10, units="cm", dpi=96)



#Cluster Dotplot of 2 imprinted genes (NOT required). 

#I took gene list names (saved as basilia_genes_for_clusterdotplot2_as_bulk.txt) to check counts
basilia_genes_for_clusterdotplot2_as_bulk <- read.table("basilia_genes_for_clusterdotplot_as_bulk2.txt", header = F)
colnames(basilia_genes_for_clusterdotplot2_as_bulk) <- c("Genes")
head(ddsNormcounts_gene)
dim(ddsNormcounts_gene)
ddsNormcounts_chr_gene1_clusterdotplot2genes <- merge(basilia_genes_for_clusterdotplot2_as_bulk,ddsNormcounts_gene, by ="Genes", all.x=FALSE,sort = F)
head(ddsNormcounts_chr_gene1_clusterdotplot2genes)
dim(ddsNormcounts_chr_gene1_clusterdotplot2genes)
library("writexl") #Remember it does not print rowname, and also I dont need rownames here in this particular  case
write_xlsx(ddsNormcounts_chr_gene1_clusterdotplot2genes,"ddsNormcounts_chr_gene1_clusterdotplot2genes.xlsx")

#First calculate average because as Basilia discussed rep1 can be compared to both rep1 and rep2 not to only one
ddsNormcounts_chr_gene1_clusterdotplot2genes["WT.avg"] <- (ddsNormcounts_chr_gene1_clusterdotplot2genes$JB1_WT_Rep1+ddsNormcounts_chr_gene1_clusterdotplot2genes$JB1_WT_Rep2)/2
ddsNormcounts_chr_gene1_clusterdotplot2genes["ZFP57KO.avg"] <- (ddsNormcounts_chr_gene1_clusterdotplot2genes$JB1_ZFP57_KO_Rep1+ddsNormcounts_chr_gene1_clusterdotplot2genes$JB1_ZFP57_KO_Rep2)/2
head(ddsNormcounts_chr_gene1_clusterdotplot2genes)

#calculate ko/(wt+ko)
ddsNormcounts_chr_gene1_clusterdotplot2genes["Ratio.avg"] <- ddsNormcounts_chr_gene1_clusterdotplot2genes$ZFP57KO.avg/(ddsNormcounts_chr_gene1_clusterdotplot2genes$WT.avg+ddsNormcounts_chr_gene1_clusterdotplot2genes$ZFP57KO.avg)
head(ddsNormcounts_chr_gene1_clusterdotplot2genes)
dim(ddsNormcounts_chr_gene1_clusterdotplot2genes)
#Separate replicates
cluster_imp_gene_indiv_bulk_part2_norm_are <- ddsNormcounts_chr_gene1_clusterdotplot2genes[,c(1,12)]
head(cluster_imp_gene_indiv_bulk_part2_norm_are)
rownames(cluster_imp_gene_indiv_bulk_part2_norm_are) <- cluster_imp_gene_indiv_bulk_part2_norm_are[,1]

cluster_imp_gene_indiv_bulk_part2_norm_sortre <- data.frame(cluster_imp_gene_indiv_bulk_part2_norm_are$Ratio.avg)
head(cluster_imp_gene_indiv_bulk_part2_norm_sortre)
rownames(cluster_imp_gene_indiv_bulk_part2_norm_sortre) <- rownames(cluster_imp_gene_indiv_bulk_part2_norm_are) 
cluster_imp_gene_indiv_bulk_part2_norm_sortre
colnames(cluster_imp_gene_indiv_bulk_part2_norm_sortre) <- "Ratio.avg"

cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot <- as.matrix(t(cluster_imp_gene_indiv_bulk_part2_norm_sortre))
head(cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot)
dim(cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot)
cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1 <- stack(cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot)
head(cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1)
colnames(cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1) <- c("Group", "Gene", "Ratio")
cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1 <- data.frame(cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1)
head(cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1)
dim(cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1)
cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1["Group.1"] <- data.frame(rep("KOWT_ratio", 2))
head(cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1)
dim(cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1)
str(cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1)
cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1_half <- cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1
cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1_half["Ratio"] <- as.numeric(cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1_half$Ratio) - 0.5
head(cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1_half)
library(ggplot2)
# Dotplot for Cluster
# Change the position
pdot2 <-ggplot(cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1_half, aes(x=Gene, y=Ratio, fill= Group.1)) + geom_hline(yintercept = c(), colour = "grey", linetype=c("dashed","solid","dashed"))+
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.8))
pdot2 + scale_fill_manual(values=c("blue"))+ylim(c(-0.5,0.5))+
  scale_color_manual(values=c("black","black","black","black")) + theme_classic()

ggsave("dotplot_cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1_half.svg", width=20, height=10, units="cm", dpi=96)
ggsave("dotplot_cluster_imp_gene_indiv_bulk_part2_norm_sortrebarplot1_half.jpg", width=20, height=10, units="cm", dpi=96)


#KS Test
#Remove duplicates due to ZFP57 peaks
sort -k5,5 -u all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrMY.txt > all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrMY_dedup.txt
sort -k5,5 -u imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrMY.txt > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrMY_dedup.txt
sort -k5,5 -u Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrMY.txt > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrMY_dedup.txt
sort -k4,4 -u res_sorted_gene_rearranged_MY.txt > res_sorted_gene_rearranged_MY_dedup.txt
sort -k4,4 -u deseq2_results_res0.05_gene_rearranged_MY.txt > deseq2_results_res0.05_gene_rearranged_MY_dedup.txt 

All_dereg_Zfp57_Kap1dedup <- read.table("all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrMY_dedup.txt")
dim(All_dereg_Zfp57_Kap1dedup) #4899 9
Imprinted_dereg_Zfp57_Kap1dedup <- read.table("imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrMY_dedup.txt")
dim(Imprinted_dereg_Zfp57_Kap1dedup) #50  9
Nonimprinted_dereg_Zfp57_Kap1dedup <- read.table("Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrMY_dedup.txt")
dim(Nonimprinted_dereg_Zfp57_Kap1dedup)#4849  9
ks.test(Imprinted_dereg_Zfp57_Kap1dedup$V9, Nonimprinted_dereg_Zfp57_Kap1dedup$V9, alternative = "greater", exact = NULL)

#For the hyper geometric test, chrM and chrY are removed. chrX was kept, deduplication was performed for all files
res_sorted_gene_rearranged_MY_dedup <- read.table("res_sorted_gene_rearranged_MY_dedup.txt")
dim(res_sorted_gene_rearranged_MY_dedup) #19955  12
deseq2_results_res0.05_gene_rearranged_MY_dedup <- read.table("deseq2_results_res0.05_gene_rearranged_MY_dedup.txt")
dim(deseq2_results_res0.05_gene_rearranged_MY_dedup) #4899   12

#Impcommand1: fgrep -f imprinted_gene_name.txt res_sorted_gene_rearranged_MY_dedup.txt -w | wc -l   123
#Impcommand2: fgrep -f imprinted_gene_name.txt deseq2_results_res0.05_gene_rearranged_MY.txt -w | wc -l 50

#Let N denote the Total number of expressed genes
#M denote the Total number of expressed imprinted genes
#K denote the Total number of DE genes
#X denote the total number of DE genes with are imprinted

#Then the p.value can be computed as

#p.value= hyper(q=X-1,m=M, n=N-M,k=K, lower.tail=FALSE)
#Hypergeomewtric Test (hgt)
hgtN = 19955 #dim(res_sorted_gene_rearranged_MY_dedup) = N, Total number of expressed genes
hgtM = 123 #see Impcommand1, ctrl+F = M, Total number of expressed imprinted genes
hgtK = 4899 #dim(deseq2_results_res0.05_gene_rearranged_MY_dedup) = K, Total number of DE genes
hgtX = 50 #see Impcommand2, ctrl+F = X,total number of DE genes with are imprinted 
hgtp.value <-  phyper(q=hgtX-1,m=hgtM, n=hgtN-hgtM,k=hgtK, lower.tail=FALSE)
hgtp.value

################### End of Analysis ###########
library(WriteXLS) or use library("writexl")
WriteXLS(ddsNormcounts,ExcelFileName="Bulk_Deseq2_Normalized_counts.xlsx", col.names=TRUE, row.names=TRUE)
WriteXLS(ddsNormcounts_gene,ExcelFileName="Bulk_Deseq2_Normalized_counts_genes.xlsx", col.names=TRUE, row.names=TRUE)
WriteXLS(ddsNormcountsavg,ExcelFileName="Bulk_Deseq2_Normalized_averagecounts_genes.xlsx", col.names=TRUE, row.names=TRUE)
WriteXLS(data.frame(res_sorted),ExcelFileName="Bulk_Deseq2_all_features.xlsx", col.names=TRUE, row.names=TRUE)
write_xlsx(deseq2_results_res0.05_gene,"Bulk_Deseq2_differentially_expressed_features.xlsx")
WriteXLS(res_sorted_gene_rearranged,ExcelFileName="Bulk_Deseq2_all_features_gene.xlsx", col.names=TRUE, row.names=TRUE)
write_xlsx(deseq2_results_res0.05_gene_rearranged,"Bulk_Deseq2_differentially_expressed_features_gene.xlsx")
WriteXLS(dgeFPKM_T12_T0_t0avg.chr.gene,ExcelFileName="Bulk_edgeR_T0_T12rep_FPKM_all_feature_gene.xlsx", col.names=TRUE, row.names=TRUE)
WriteXLS(selected_markers_T0_T12_1.norm,ExcelFileName="Marker_analysis_25markers.xlsx", col.names=TRUE, row.names=TRUE)
WriteXLS(scatterund_bulk_mm10_t12_column,ExcelFileName="Scatter_plot_log_t12.xlsx", col.names=TRUE, row.names=TRUE)
WriteXLS(gost.diff_0.5.up.res.sorted,ExcelFileName="Bulk_Gprofiler_Upregulated_Genes.xlsx", col.names=TRUE, row.names=TRUE)
WriteXLS(gost.diff_0.5.down.res.sorted,ExcelFileName="Bulk_Gprofiler_Downregulated_Genes.xlsx", col.names=TRUE, row.names=TRUE)
write_xlsx(gost.diff_0.5.downun100kb.res.sorted,"Bulk_Gprofiler_Downregulatedunder100kb_Genes.xlsx")
WriteXLS(plotdistance,ExcelFileName="Bulk_distanceplot_values.xlsx", col.names=TRUE, row.names=TRUE)
write_xlsx(scatterund_bulk_mm10_t12_column_uniq,"Scatter_plot_uniq_log_t12.xlsx")


#All Possible type Functional analysis
#Extract deregulated genes underlying Zfp57 peaks windows
sort -k5,5 imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -u | awk '{print $5}' > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.gene
sort -k5,5 Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -u | awk '{print $5}' > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.gene

diff.imprinted.symbol <- read.table("imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.gene", header =F)
diff.imprinted.symbol <-diff.imprinted.symbol$V1
gost.diff.imprinted <- gost(query=as.character(diff.imprinted.symbol), organism="mmusculus")
gostplot(gost.diff.imprinted)
#Sort by p-value
gost.diff.imprinted.res.sorted <- gost.diff.imprinted$result[order(gost.diff.imprinted$result$p_value),]
head(gost.diff.imprinted.res.sorted)
gost.diff.imprinted.res.sorted <- gost.diff.imprinted$result[order(gost.diff.imprinted$result$p_value),]
head(gost.diff.imprinted.res.sorted)

#GO:BP
gost.diff.imprinted.res.sorted_gobp <- gost.diff.imprinted.res.sorted[which(gost.diff.imprinted.res.sorted$source == "GO:BP"),]
head(gost.diff.imprinted.res.sorted_gobp)
dim(gost.diff.imprinted.res.sorted_gobp)
gost.diff.imprinted.res.sorted_gobp_bar <- gost.diff.imprinted.res.sorted_gobp[,c(11,3)]
head(gost.diff.imprinted.res.sorted_gobp_bar)
gost.diff.imprinted.res.sorted_gobp_bar_top <- head(gost.diff.imprinted.res.sorted_gobp_bar,10)
gost.diff.imprinted.res.sorted_gobp_bar_top$term_name <- gsub(' ', '.', gost.diff.imprinted.res.sorted_gobp_bar_top$term_name)
gost.diff.imprinted.res.sorted_gobp_bar_top$p_value <- -log10(gost.diff.imprinted.res.sorted_gobp_bar_top$p_value)
gost.diff.imprinted.res.sorted_gobp_bar_top <- data.frame(gost.diff.imprinted.res.sorted_gobp_bar_top)

library(ggpubr)
ggbarplot(gost.diff.imprinted.res.sorted_gobp_bar_top, 
          x = "term_name", 
          y = "p_value", 
          color = "#e58e8eff", 
          fill = "#e58e8eff" ,
          sort.by.grodowns = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "-log10(p.value)",
          xlab = "Pathways",
          legend.title = "Pathways",
          lab.size = 9,
          rotate = TRUE, position = position_dodge(),
          sort.val = "asc",
          ggtheme = theme_bw())

ggsave("gost.diff.imprinted.res.sorted_gobp_bar_top.svg", width=13.5, height=8, units="cm", dpi=96)


diff.Nonimprinted.symbol <- read.table("Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.gene", header =F)
diff.Nonimprinted.symbol <-diff.Nonimprinted.symbol$V1
gost.diff.Nonimprinted <- gost(query=as.character(diff.Nonimprinted.symbol), organism="mmusculus")
gostplot(gost.diff.Nonimprinted)
#Sort by p-value
gost.diff.Nonimprinted.res.sorted <- gost.diff.Nonimprinted$result[order(gost.diff.Nonimprinted$result$p_value),]
head(gost.diff.Nonimprinted.res.sorted)
gost.diff.Nonimprinted.res.sorted <- gost.diff.Nonimprinted$result[order(gost.diff.Nonimprinted$result$p_value),]
head(gost.diff.Nonimprinted.res.sorted)

#GO:BP
gost.diff.Nonimprinted.res.sorted_gobp <- gost.diff.Nonimprinted.res.sorted[which(gost.diff.Nonimprinted.res.sorted$source == "GO:BP"),]
head(gost.diff.Nonimprinted.res.sorted_gobp)
dim(gost.diff.Nonimprinted.res.sorted_gobp)
gost.diff.Nonimprinted.res.sorted_gobp_bar <- gost.diff.Nonimprinted.res.sorted_gobp[,c(11,3)]
head(gost.diff.Nonimprinted.res.sorted_gobp_bar)
gost.diff.Nonimprinted.res.sorted_gobp_bar_top <- head(gost.diff.Nonimprinted.res.sorted_gobp_bar,10)
gost.diff.Nonimprinted.res.sorted_gobp_bar_top$term_name <- gsub(' ', '.', gost.diff.Nonimprinted.res.sorted_gobp_bar_top$term_name)
gost.diff.Nonimprinted.res.sorted_gobp_bar_top$p_value <- -log10(gost.diff.Nonimprinted.res.sorted_gobp_bar_top$p_value)
gost.diff.Nonimprinted.res.sorted_gobp_bar_top <- data.frame(gost.diff.Nonimprinted.res.sorted_gobp_bar_top)

library(ggpubr)
ggbarplot(gost.diff.Nonimprinted.res.sorted_gobp_bar_top, 
          x = "term_name", 
          y = "p_value", 
          color = "#e58e8eff", 
          fill = "#e58e8eff" ,
          sort.by.grodowns = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "-log10(p.value)",
          xlab = "Pathways",
          legend.title = "Pathways",
          lab.size = 9,
          rotate = TRUE, position = position_dodge(),
          sort.val = "asc",
          ggtheme = theme_bw())

ggsave("gost.diff.Nonimprinted.res.sorted_gobp_bar_top.svg", width=13.5, height=8, units="cm", dpi=96)


awk '$9 <= 100000 {print $0}' imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt | sort -k5,5 -u > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_100000.txt
awk '$9 <= 100000 {print $0}' Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt | sort -k5,5 -u > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_100000.txt

#Extract deregulated genes underlying Zfp57 peaks windows
sort -k5,5 imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_100000.txt -u | awk '{print $5}' > imprinted100Kb_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.gene
sort -k5,5 Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_100000.txt -u | awk '{print $5}' > Nonimprinted100Kb_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.gene

diff.imprinted100Kb.symbol <- read.table("imprinted100Kb_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.gene", header =F)
diff.imprinted100Kb.symbol <-diff.imprinted100Kb.symbol$V1
gost.diff.imprinted100Kb <- gost(query=as.character(diff.imprinted100Kb.symbol), organism="mmusculus")
gostplot(gost.diff.imprinted100Kb)
#Sort by p-value
gost.diff.imprinted100Kb.res.sorted <- gost.diff.imprinted100Kb$result[order(gost.diff.imprinted100Kb$result$p_value),]
head(gost.diff.imprinted100Kb.res.sorted)
gost.diff.imprinted100Kb.res.sorted <- gost.diff.imprinted100Kb$result[order(gost.diff.imprinted100Kb$result$p_value),]
head(gost.diff.imprinted100Kb.res.sorted)

#GO:BP
gost.diff.imprinted100Kb.res.sorted_gobp <- gost.diff.imprinted100Kb.res.sorted[which(gost.diff.imprinted100Kb.res.sorted$source == "GO:BP"),]
head(gost.diff.imprinted100Kb.res.sorted_gobp)
dim(gost.diff.imprinted100Kb.res.sorted_gobp)
gost.diff.imprinted100Kb.res.sorted_gobp_bar <- gost.diff.imprinted100Kb.res.sorted_gobp[,c(11,3)]
head(gost.diff.imprinted100Kb.res.sorted_gobp_bar)
gost.diff.imprinted100Kb.res.sorted_gobp_bar_top <- head(gost.diff.imprinted100Kb.res.sorted_gobp_bar,10)
gost.diff.imprinted100Kb.res.sorted_gobp_bar_top$term_name <- gsub(' ', '.', gost.diff.imprinted100Kb.res.sorted_gobp_bar_top$term_name)
gost.diff.imprinted100Kb.res.sorted_gobp_bar_top$p_value <- -log10(gost.diff.imprinted100Kb.res.sorted_gobp_bar_top$p_value)
gost.diff.imprinted100Kb.res.sorted_gobp_bar_top <- data.frame(gost.diff.imprinted100Kb.res.sorted_gobp_bar_top)

library(ggpubr)
ggbarplot(gost.diff.imprinted100Kb.res.sorted_gobp_bar_top, 
          x = "term_name", 
          y = "p_value", 
          color = "#e58e8eff", 
          fill = "#e58e8eff" ,
          sort.by.grodowns = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "-log10(p.value)",
          xlab = "Pathways",
          legend.title = "Pathways",
          lab.size = 9,
          rotate = TRUE, position = position_dodge(),
          sort.val = "asc",
          ggtheme = theme_bw())

ggsave("gost.diff.imprinted100Kb.res.sorted_gobp_bar_top.svg", width=13.5, height=8, units="cm", dpi=96)


diff.Nonimprinted100Kb.symbol <- read.table("Nonimprinted100Kb_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.gene", header =F)
diff.Nonimprinted100Kb.symbol <-diff.Nonimprinted100Kb.symbol$V1
gost.diff.Nonimprinted100Kb <- gost(query=as.character(diff.Nonimprinted100Kb.symbol), organism="mmusculus")
gostplot(gost.diff.Nonimprinted100Kb)
#Sort by p-value
gost.diff.Nonimprinted100Kb.res.sorted <- gost.diff.Nonimprinted100Kb$result[order(gost.diff.Nonimprinted100Kb$result$p_value),]
head(gost.diff.Nonimprinted100Kb.res.sorted)
gost.diff.Nonimprinted100Kb.res.sorted <- gost.diff.Nonimprinted100Kb$result[order(gost.diff.Nonimprinted100Kb$result$p_value),]
head(gost.diff.Nonimprinted100Kb.res.sorted)

#GO:BP
gost.diff.Nonimprinted100Kb.res.sorted_gobp <- gost.diff.Nonimprinted100Kb.res.sorted[which(gost.diff.Nonimprinted100Kb.res.sorted$source == "GO:BP"),]
head(gost.diff.Nonimprinted100Kb.res.sorted_gobp)
dim(gost.diff.Nonimprinted100Kb.res.sorted_gobp)
gost.diff.Nonimprinted100Kb.res.sorted_gobp_bar <- gost.diff.Nonimprinted100Kb.res.sorted_gobp[,c(11,3)]
head(gost.diff.Nonimprinted100Kb.res.sorted_gobp_bar)
gost.diff.Nonimprinted100Kb.res.sorted_gobp_bar_top <- head(gost.diff.Nonimprinted100Kb.res.sorted_gobp_bar,10)
gost.diff.Nonimprinted100Kb.res.sorted_gobp_bar_top$term_name <- gsub(' ', '.', gost.diff.Nonimprinted100Kb.res.sorted_gobp_bar_top$term_name)
gost.diff.Nonimprinted100Kb.res.sorted_gobp_bar_top$p_value <- -log10(gost.diff.Nonimprinted100Kb.res.sorted_gobp_bar_top$p_value)
gost.diff.Nonimprinted100Kb.res.sorted_gobp_bar_top <- data.frame(gost.diff.Nonimprinted100Kb.res.sorted_gobp_bar_top)

library(ggpubr)
ggbarplot(gost.diff.Nonimprinted100Kb.res.sorted_gobp_bar_top, 
          x = "term_name", 
          y = "p_value", 
          color = "#e58e8eff", 
          fill = "#e58e8eff" ,
          sort.by.grodowns = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "-log10(p.value)",
          xlab = "Pathways",
          legend.title = "Pathways",
          lab.size = 9,
          rotate = TRUE, position = position_dodge(),
          sort.val = "asc",
          ggtheme = theme_bw())

ggsave("gost.diff.Nonimprinted100Kb.res.sorted_gobp_bar_top.svg", width=13.5, height=8, units="cm", dpi=96)

awk '$9 <= 1000000 {print $0}' imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt | sort -k5,5 -u > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_1000000.txt
awk '$9 <= 1000000 {print $0}' Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt | sort -k5,5 -u > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_1000000.txt

#Extract deregulated genes underlying Zfp57 peaks windows
sort -k5,5 imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_1000000.txt -u | awk '{print $5}' > imprinted1Mb_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.gene
sort -k5,5 Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_1000000.txt -u | awk '{print $5}' > Nonimprinted1Mb_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.gene

diff.imprinted1Mb.symbol <- read.table("imprinted1Mb_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.gene", header =F)
diff.imprinted1Mb.symbol <-diff.imprinted1Mb.symbol$V1
gost.diff.imprinted1Mb <- gost(query=as.character(diff.imprinted1Mb.symbol), organism="mmusculus")
gostplot(gost.diff.imprinted1Mb)
#Sort by p-value
gost.diff.imprinted1Mb.res.sorted <- gost.diff.imprinted1Mb$result[order(gost.diff.imprinted1Mb$result$p_value),]
head(gost.diff.imprinted1Mb.res.sorted)
gost.diff.imprinted1Mb.res.sorted <- gost.diff.imprinted1Mb$result[order(gost.diff.imprinted1Mb$result$p_value),]
head(gost.diff.imprinted1Mb.res.sorted)

#GO:BP
gost.diff.imprinted1Mb.res.sorted_gobp <- gost.diff.imprinted1Mb.res.sorted[which(gost.diff.imprinted1Mb.res.sorted$source == "GO:BP"),]
head(gost.diff.imprinted1Mb.res.sorted_gobp)
dim(gost.diff.imprinted1Mb.res.sorted_gobp)
gost.diff.imprinted1Mb.res.sorted_gobp_bar <- gost.diff.imprinted1Mb.res.sorted_gobp[,c(11,3)]
head(gost.diff.imprinted1Mb.res.sorted_gobp_bar)
gost.diff.imprinted1Mb.res.sorted_gobp_bar_top <- head(gost.diff.imprinted1Mb.res.sorted_gobp_bar,10)
gost.diff.imprinted1Mb.res.sorted_gobp_bar_top$term_name <- gsub(' ', '.', gost.diff.imprinted1Mb.res.sorted_gobp_bar_top$term_name)
gost.diff.imprinted1Mb.res.sorted_gobp_bar_top$p_value <- -log10(gost.diff.imprinted1Mb.res.sorted_gobp_bar_top$p_value)
gost.diff.imprinted1Mb.res.sorted_gobp_bar_top <- data.frame(gost.diff.imprinted1Mb.res.sorted_gobp_bar_top)

library(ggpubr)
ggbarplot(gost.diff.imprinted1Mb.res.sorted_gobp_bar_top, 
          x = "term_name", 
          y = "p_value", 
          color = "#e58e8eff", 
          fill = "#e58e8eff" ,
          sort.by.grodowns = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "-log10(p.value)",
          xlab = "Pathways",
          legend.title = "Pathways",
          lab.size = 9,
          rotate = TRUE, position = position_dodge(),
          sort.val = "asc",
          ggtheme = theme_bw())

ggsave("gost.diff.imprinted1Mb.res.sorted_gobp_bar_top.svg", width=13.5, height=8, units="cm", dpi=96)


diff.Nonimprinted1Mb.symbol <- read.table("Nonimprinted1Mb_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.gene", header =F)
diff.Nonimprinted1Mb.symbol <-diff.Nonimprinted1Mb.symbol$V1
gost.diff.Nonimprinted1Mb <- gost(query=as.character(diff.Nonimprinted1Mb.symbol), organism="mmusculus")
gostplot(gost.diff.Nonimprinted1Mb)
#Sort by p-value
gost.diff.Nonimprinted1Mb.res.sorted <- gost.diff.Nonimprinted1Mb$result[order(gost.diff.Nonimprinted1Mb$result$p_value),]
head(gost.diff.Nonimprinted1Mb.res.sorted)
gost.diff.Nonimprinted1Mb.res.sorted <- gost.diff.Nonimprinted1Mb$result[order(gost.diff.Nonimprinted1Mb$result$p_value),]
head(gost.diff.Nonimprinted1Mb.res.sorted)

#GO:BP
gost.diff.Nonimprinted1Mb.res.sorted_gobp <- gost.diff.Nonimprinted1Mb.res.sorted[which(gost.diff.Nonimprinted1Mb.res.sorted$source == "GO:BP"),]
head(gost.diff.Nonimprinted1Mb.res.sorted_gobp)
dim(gost.diff.Nonimprinted1Mb.res.sorted_gobp)
gost.diff.Nonimprinted1Mb.res.sorted_gobp_bar <- gost.diff.Nonimprinted1Mb.res.sorted_gobp[,c(11,3)]
head(gost.diff.Nonimprinted1Mb.res.sorted_gobp_bar)
gost.diff.Nonimprinted1Mb.res.sorted_gobp_bar_top <- head(gost.diff.Nonimprinted1Mb.res.sorted_gobp_bar,10)
gost.diff.Nonimprinted1Mb.res.sorted_gobp_bar_top$term_name <- gsub(' ', '.', gost.diff.Nonimprinted1Mb.res.sorted_gobp_bar_top$term_name)
gost.diff.Nonimprinted1Mb.res.sorted_gobp_bar_top$p_value <- -log10(gost.diff.Nonimprinted1Mb.res.sorted_gobp_bar_top$p_value)
gost.diff.Nonimprinted1Mb.res.sorted_gobp_bar_top <- data.frame(gost.diff.Nonimprinted1Mb.res.sorted_gobp_bar_top)

library(ggpubr)
ggbarplot(gost.diff.Nonimprinted1Mb.res.sorted_gobp_bar_top, 
          x = "term_name", 
          y = "p_value", 
          color = "#e58e8eff", 
          fill = "#e58e8eff" ,
          sort.by.grodowns = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "-log10(p.value)",
          xlab = "Pathways",
          legend.title = "Pathways",
          lab.size = 9,
          rotate = TRUE, position = position_dodge(),
          sort.val = "asc",
          ggtheme = theme_bw())

ggsave("gost.diff.Nonimprinted1Mb.res.sorted_gobp_bar_top.svg", width=13.5, height=8, units="cm", dpi=96)









######################################## extra commands ################
#Plot PCA
plotPCA(as.matrix(ddsNormcounts))  #library(EDASeq)
plotPCA(as.matrix(ddsNormcounts), labels=F, col =  c("blue","cyan","darkred","magenta"))

ggsave("PCA_tddsNormcountsfilt_scaleT.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)


tddsNormcountsfilt = t(ddsNormcounts)
dim(tddsNormcountsfilt)
rownames(tddsNormcountsfilt)
tddsNormcountsfilt = data.frame(tddsNormcountsfilt)
tddsNormcountsfilt["Color"] <-  c("WT_1","WT_2","ZFP57KO_1","ZFP57KO_2")
dim(tddsNormcountsfilt)
dfx <-tddsNormcountsfilt[c(1:19995)]
PC<-prcomp(dfx, scale. = T)
PCi<-data.frame(PC$x,Color=tddsNormcountsfilt$Color)
percentage <- round(PC$sdev^2 / sum(PC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentage <- paste( colnames(PCi), "(", paste( as.character(percentage), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

p15<-ggplot(PCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentage[1]) + ylab(percentage[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("blue","cyan","darkred","magenta"))+
  scale_shape_manual(values=c(19,17,19,17,19,17,19,17))
p15 <- p15+theme_bw()
p15

ggsave("PCA_tddsNormcountsfilt_scaleT_man.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

############################################## NOISeq ####################################################àà

library(NOISeq)

setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2")
countdata <- read.table("JB1_time12_star-featureCounts_mm10.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
#Extract length: 
mylength <- countdata[ ,5]
countdata <- countdata[ ,6:ncol(countdata)]
head(countdata)
# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("\\.[sb]am$", "", colnames(countdata))

# Convert to matrix
mycounts <- as.matrix(countdata)
head(mycounts)
rownames(mycounts)
colnames(mycounts)
dim(mycounts)
myfactors = read.table("myfactors.txt", header=TRUE)
head(myfactors)
rownames(myfactors)
colnames(myfactors)

#Plot first PCA without filtering or normalisation of data: see ggplot2 
mydata1 <- NOISeq::readData(data=mycounts, factors=myfactors)
myPCA1 = dat(mydata1, type = "PCA")
par(mfrow = c(1, 2))
explo.plot(myPCA1, factor = "Sample")

#Low-count filtering
#You could use filtered.data first to remove low count features across all samples, and then use noiseqbio with the argument filter = 0 so that it does not perform any filtering.
myfilt = filtered.data(mycounts, factor = myfactors$Sample, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1, p.adj = "fdr")
head(myfilt)
mydata2 <- NOISeq::readData(data=myfilt, factors=myfactors)
myPCA2 = dat(mydata2, type = "PCA")
par(mfrow = c(1, 2))
explo.plot(myPCA2, factor = "Sample")

#Normalisation
myTMM = tmm(assayData(mydata2)$exprs, long = 1000, lc = 0)
head(myTMM)
mydata3 <- NOISeq::readData(data=myTMM, factors=myfactors)
myPCA3 = dat(mydata3, type = "PCA")
par(mfrow = c(1, 2))
explo.plot(myPCA3, factor = "Sample")


#Quality Control report
QCreport(mydata3, samples = NULL, factor = "Sample", norm = FALSE)

mydatax = mydata3
dim(mydata3)

#Sort mycounts
mycounts_sorted = mycounts[order(rownames(mycounts)),]
head(mycounts_sorted)
mycounts_sorted_sum = cbind(mycounts_sorted, mycounts_sorted[,1] + mycounts_sorted[,2],
                            mycounts_sorted[,3] + mycounts_sorted[,4])


head(mycounts_sorted_sum)
dim(mycounts_sorted_sum)

#Validation
head(mycounts)
write.table(mycounts, "mycounts.txt", sep="\t", quote = FALSE, append = FALSE)
#Validation
head(myTMM)
write.table(myTMM, "mytmm.txt", sep="\t", quote = FALSE, append = FALSE)
#Validation
head(myfilt)
write.table(myfilt, "myfilt.txt", sep="\t", quote = FALSE, append = FALSE)

mynoiseqbio = noiseqbio(mydatax, k = 0.5, norm = "n", factor = "Sample", conditions = c("WT","KO"), lc = 0, r = 50, adj = 1.5, plot = TRUE, a0per = 0.9, random.seed = 12345, filter = 0)
deg = degenes(mynoiseqbio, q = 0.98, M = NULL)
deg_sorted = deg[order(rownames(deg)),]
head(deg_sorted)
write.table(deg_sorted, "mynoiseqbio.deg_sorted.txt", sep="\t", quote = FALSE, append = FALSE)
DE.plot(mynoiseqbio, q = 0.98, graphic = "expr", log.scale = TRUE)
DE.plot(mynoiseqbio, q = 0.98, graphic = "MD")


mynoiseqbiorpkm = noiseqbio(mydata1, k = 0.5, norm = "rpkm", factor="Sample", lc = 0, r = 50, adj = 1.5, plot = FALSE,   a0per = 0.9, random.seed = 12345, filter = 1)
degrpkm = degenes(mynoiseqbiorpkm, q = 0.98, M = NULL)
degrpkm_sorted = degrpkm[order(rownames(degrpkm)),]
head(degrpkm_sorted)
write.table(degrpkm_sorted, "mynoiseqbio.degrpkm_sorted.txt", sep="\t", quote = FALSE, append = FALSE)



#------------Overlapped Strategy-----------------------------------------------------------------------------------------------
#Import ensembl ID of deseq2_mynoiseq_overlapped_features
deseq2_mynoiseq_overlapped_0.05_features <- read.table("deseq2_mynoiseq_overlapped_0.05_features.txt" , header = FALSE, stringsAsFactors = FALSE)
head(deseq2_mynoiseq_overlapped_0.05_features)
colnames(deseq2_mynoiseq_overlapped_0.05_features) <- "id"
colnames(deseq2_mynoiseq_overlapped_0.05_features)
dim(deseq2_mynoiseq_overlapped_0.05_features)
deseq2_mynoiseq_overlapped_features <- read.table("deseq2_mynoiseq_overlapped_features.txt" , header = FALSE, stringsAsFactors = FALSE)
head(deseq2_mynoiseq_overlapped_features)
colnames(deseq2_mynoiseq_overlapped_features) <- "id"
colnames(deseq2_mynoiseq_overlapped_features)
dim(deseq2_mynoiseq_overlapped_features)

#Noiseq outputs
#deg_WT_KO_sorted
#Deseq2 outputs
#res_WT_KO_sorted
#0.05
#deseq2_results_res_WT_KO_0.05

#----------------------------------------Noiseq outputs------------------------------------------------------------
#deg_sorted
head(deg_sorted)
deg_sorted1 <- cbind(rownames(deg_sorted),deg_sorted[,1:5])
head(deg_sorted1)
colnames(deg_sorted1) <- c("id","KO_mean", "WT_mean","theta","prob","log2FC")
head(deg_sorted1)
deg_sorted1 = merge(deseq2_mynoiseq_overlapped_features, deg_sorted1, by="id", all.x=TRUE)
head(deg_sorted1)
deg_sorted1_0.05 = merge(deseq2_mynoiseq_overlapped_0.05_features, deg_sorted1, by="id", all.x=TRUE)
head(deg_sorted1_0.05)

mynoiseqbio.merge.statistics <- cbind(deg_sorted1)
head(mynoiseqbio.merge.statistics)
write.table(mynoiseqbio.merge.statistics, "mynoiseqbio.merge.statistics.txt", sep="\t", quote = FALSE, append = FALSE)

mynoiseqbio.merge.statistics_0.05 <- cbind(deg_sorted1_0.05)
head(mynoiseqbio.merge.statistics_0.05)
write.table(mynoiseqbio.merge.statistics_0.05, "mynoiseqbio.merge.statistics.0.05.txt", sep="\t", quote = FALSE, append = FALSE)

#Deseq2 outputs-------------No Deseq2 pdj cutoffs---------------------------------------------------------------------------------------
#res_sorted
head(res_sorted)
res_sorted_a <- cbind(rownames(res_sorted),res_sorted[,1:7])
head(res_sorted_a)
colnames(res_sorted_a) <- c("id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue","padj", "threshold")
head(res_sorted_a)
res_sorted_a = merge(deseq2_mynoiseq_overlapped_features, res_sorted_a, by="id", all.x=TRUE)
head(res_sorted_a)

deseq2.merge.statistics <- cbind(res_sorted_a)
head(deseq2.merge.statistics)
dim(deseq2.merge.statistics)
write.table(deseq2.merge.statistics, "deseq2.merge.statistics.txt", sep="\t", quote = FALSE, append = FALSE)

#Deseq2 outputs-------------Deseq2 pdj 0.05 cutoffs---------------------------------------------------------------------------------------
#deseq2_results_res_0.05
head(deseq2_results_res0.05)
deseq2_results_res0.05_a <- cbind(rownames(deseq2_results_res0.05),deseq2_results_res0.05[,1:7])
head(deseq2_results_res0.05_a)
colnames(deseq2_results_res0.05_a) <- c("id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue","padj", "threshold")
head(deseq2_results_res0.05_a)
deseq2_results_res0.05_a = merge(deseq2_mynoiseq_overlapped_0.05_features, deseq2_results_res0.05_a, by="id", all.x=TRUE)
head(deseq2_results_res0.05_a)
dim(deseq2_results_res0.05_a)
deseq2.merge.0.05.statistics <- cbind(deseq2_results_res0.05_a)
head(deseq2.merge.0.05.statistics)
write.table(deseq2.merge.0.05.statistics, "deseq2.merge.0.05.statistics.txt", sep="\t", quote = FALSE, append = FALSE)


#################### NON OVERLAPPED STRATEGY: Only DESeq2 outputs: ALL filtered gene, no-cutoff ##########################
#Chromosome positions
chr.pos =  read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt",header=FALSE)
head(chr.pos)
colnames(chr.pos) <- c("id", "Genes", "chr", "start", "end")
head(chr.pos)
chr.pos.1 = merge(deseq2_mynoiseq_overlapped_features, chr.pos, by="id", all.x=TRUE)
head(chr.pos.1)
dim(chr.pos.1)
head(mycounts_sorted_sum)
mycounts_sorted_sum_id <- data.frame(cbind(rownames(mycounts_sorted_sum), mycounts_sorted_sum))
head(mycounts_sorted_sum_id)


colnames(mycounts_sorted_sum_id) <- c("id","JB1_WT_Rep1","JB1_WT_Rep2","JB1_ZFP57_KO_Rep1","JB1_ZFP57_KO_Rep2","JB1_WT_Rep.CPM","JB1_ZFP57_KO_Rep.CPM")
head(mycounts_sorted_sum_id)
dim(mycounts_sorted_sum_id)
mycounts_sorted_sum_id_re = merge(deseq2_mynoiseq_overlapped_features, mycounts_sorted_sum_id, by="id", all.x=TRUE)
head(mycounts_sorted_sum_id_re)
dim(mycounts_sorted_sum_id_re)
mycounts_sorted_sum_genechr <- data.frame(cbind(chr.pos.1, mycounts_sorted_sum_id_re))
head(mycounts_sorted_sum_genechr)
mycounts_noiseeq_deseq2_sorted_sum_genechr <- data.frame(cbind(chr.pos.1, mycounts_sorted_sum_id_re, mynoiseqbio.merge.statistics, deseq2.merge.statistics))
head(mycounts_noiseeq_deseq2_sorted_sum_genechr)
write.table(mycounts_noiseeq_deseq2_sorted_sum_genechr, "mycounts_noiseeq_deseq2_sorted_sum_genechr.txt", sep="\t", quote = FALSE, append = FALSE)

#Extract mycounts counts and Prepare CPM normalizations
mycounts_sorted_sum_col <- mycounts_sorted_sum
colnames(mycounts_sorted_sum_col) <- c("JB1_WT_Rep1","JB1_WT_Rep2","JB1_ZFP57_KO_Rep1","JB1_ZFP57_KO_Rep2","JB1_WT_Rep.CPM","JB1_ZFP57_KO_Rep.CPM")
head(mycounts_sorted_sum_col)
dim(mycounts_sorted_sum_col)
colSums(mycounts_sorted_sum_col)
mycounts_sorted_CPM = cbind(mycounts_sorted, ((mycounts_sorted[,1] + mycounts_sorted[,2])*1000000)/17884631,
                            ((mycounts_sorted[,3] + mycounts_sorted[,4])*1000000)/19881558)


head(mycounts_sorted_CPM)
mycounts_sorted_CPM_id <- data.frame(cbind(rownames(mycounts_sorted_CPM), mycounts_sorted_CPM))
head(mycounts_sorted_CPM_id)
colnames(mycounts_sorted_CPM_id) <- c("id","JB1_WT_Rep1","JB1_WT_Rep2","JB1_ZFP57_KO_Rep1","JB1_ZFP57_KO_Rep2","JB1_WT_Rep.CPM","JB1_ZFP57_KO_Rep.CPM")
head(mycounts_sorted_CPM_id)
dim(mycounts_sorted_CPM_id)
mycounts_sorted_CPM_id_re = merge(deseq2_mynoiseq_overlapped_features, mycounts_sorted_CPM_id, by="id", all.x=TRUE)
head(mycounts_sorted_CPM_id_re)
mycounts_sorted_CPM_genechr <- data.frame(cbind(chr.pos.1, mycounts_sorted_CPM_id_re))
head(mycounts_sorted_CPM_genechr)
mycounts_noiseeq_deseq2_sorted_CPM_genechr <- data.frame(cbind(chr.pos.1, mycounts_sorted_CPM_id_re, mynoiseqbio.merge.statistics, deseq2.merge.statistics))
head(mycounts_noiseeq_deseq2_sorted_CPM_genechr,2)
write.table(mycounts_noiseeq_deseq2_sorted_CPM_genechr, "mycounts_noiseeq_deseq2_sorted_CPM_genechr.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)
mycounts_sorted_CPM_genechr_rearranged <- mycounts_sorted_CPM_genechr[,c(3:5, 2,1,6:12)]
head(mycounts_sorted_CPM_genechr_rearranged,2)
mycounts_noiseeq_deseq2_sorted_CPM_genechr_rearranged <- data.frame(cbind(mycounts_sorted_CPM_genechr_rearranged, mynoiseqbio.merge.statistics, deseq2.merge.statistics))
head(mycounts_noiseeq_deseq2_sorted_CPM_genechr_rearranged,2)
write.table(mycounts_noiseeq_deseq2_sorted_CPM_genechr_rearranged, "mycounts_noiseeq_deseq2_sorted_CPM_genechr_rearranged.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)


#z-score calc
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/Bulk_mm10/featurecount")
ddsNormcounts <- read.table("ddsNormcounts.txt" , header = TRUE, stringsAsFactors = FALSE)
head(ddsNormcounts)
dim(ddsNormcounts)
colnames(ddsNormcounts) <- c("JB1_WT_Rep1","JB1_WT_Rep2","JB1_ZFP57_KO_Rep1","JB1_ZFP57_KO_Rep2")
ddsNormcounts <- data.frame(ddsNormcounts[order(rownames(ddsNormcounts)),])
head(ddsNormcounts)

#Transform and scale
z_TddsNormcounts= scale(t(ddsNormcounts), center = TRUE, scale = TRUE)
z_ddsNormcounts <- t(z_TddsNormcounts)
write.table(z_ddsNormcounts, "z_ddsNormcounts.txt", sep="\t", quote=F, col.names=NA)
head(z_ddsNormcounts)
#z-score calc on averaged normalized counts
head(ddsNormcounts)
dim(ddsNormcounts)
#ddsNormcountsavg <- data.frame(cbind(rownames(ddsNormcounts),(ddsNormcounts$JB1_WT_Rep1+ddsNormcounts$JB1_WT_Rep2)/2, (ddsNormcounts$JB1_ZFP57_KO_Rep1 + ddsNormcounts$JB1_ZFP57_KO_Rep2)/2))
ddsNormcountsavg <- data.frame(cbind((ddsNormcounts$JB1_WT_Rep1+ddsNormcounts$JB1_WT_Rep2)/2, (ddsNormcounts$JB1_ZFP57_KO_Rep1 + ddsNormcounts$JB1_ZFP57_KO_Rep2)/2))

rownames(ddsNormcountsavg) <- rownames(ddsNormcounts)
head(ddsNormcountsavg)
dim(ddsNormcountsavg)
colnames(ddsNormcountsavg) <- c("JB1_WT","JB1_ZFP57KO")
rownames(ddsNormcountsavg)
head(ddsNormcountsavg)

#Transform and scale
z_TddsNormcountsavg= scale(t(ddsNormcountsavg), center = TRUE, scale = TRUE)
z_ddsNormcountsavg <- t(z_TddsNormcountsavg)
head(z_ddsNormcountsavg)
dim(z_ddsNormcountsavg)
write.table(z_ddsNormcountsavg, "z_ddsNormcountsavg.txt", sep="\t", quote=F, col.names=NA)

#Without trnasformation:column wise z-scoring
Nz_TddsNormcountsavg = scale((ddsNormcountsavg), center = TRUE, scale = TRUE)
Nz_ddsNormcountsavg <- Nz_TddsNormcountsavg
write.table(Nz_ddsNormcountsavg, "Nz_ddsNormcountsavg.txt", sep="\t", quote=F, col.names=NA)
head(Nz_ddsNormcountsavg)


#Without trnasformation:column wise z-scoring
NJB1_T12_fpkm_matrixavg = scale((JB1_T12_fpkm_matrixavg), center = TRUE, scale = TRUE)
head(NJB1_T12_fpkm_matrixavg)
##write.table(NJB1_T12_fpkm_matrixavg, "NJB1_T12_fpkm_matrixavg.txt", sep="\t", quote=F, col.names=NA)
##write.table(NJB1_T12_fpkm_matrixavg, "NJB1_T12_fpkm_matrix_GeneLengthE_S_zscore.txt", sep="\t", quote = FALSE, append = FALSE, row.names = TRUE)
write.table(NJB1_T12_fpkm_matrixavg, "NJB1_T12_fpkm_matrix_GeneLengthFeatureCount_zscore.txt", sep="\t", quote = FALSE, append = FALSE, row.names = TRUE)



#Calculate RPKM/FPKM by NOIseq
library(NOISeq)
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2")
countdata <- read.table("JB1_time12_star-featureCounts_mm10.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
#Extract length: 
mylength <- read.table("mylength.txt", header=TRUE)
head(mylength)

countdata <- countdata[ ,6:ncol(countdata)]
head(countdata)
# Remove .bam or .sam from filenames
colnames(countdata) <- gsub("_Aligned.sortedByCoord.out.bam", "", colnames(countdata))
head(countdata)
dim(countdata)

# Convert to matrix
mycounts <- as.matrix(countdata)
head(mycounts)
rownames(mycounts)
colnames(mycounts)
dim(mycounts)


#################################################### Transform data #####################################à
head(ddscounts)
ddscounts1 <- as.matrix(ddscounts)
tddscounts1 = t(ddscounts1)
dim(tddscounts1)
rownames(tddscounts1)
tddscounts1 = data.frame(tddscounts1)
write.table(tddscounts1 , "tddscounts1.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tddscounts1)
tddscounts1["Color"] <-  c("wt","wt","ko","ko")
dfx <-tddscounts1[c(1:19995)]
PC<-prcomp(dfx, scale. = T)
head(PC)
PCi<-data.frame(PC$x,Color=tddscounts1$Color)
percentage <- round(PC$sdev^2 / sum(PC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentage <- paste( colnames(PCi), "(", paste( as.character(percentage), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

p3<-ggplot(PCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentage[1]) + ylab(percentage[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("darkred","darkgreen"))+
  scale_shape_manual(values=c(19,19,19,19,19,19))
p3 <- p3+theme_classic()
p3
ggsave("PCA_ddscounts.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96)


#################################################### Transform data #####################################à
head(ddsNormcounts)
ddsNormcounts1 <- as.matrix(ddsNormcounts)
tddsNormcounts1 = t(ddsNormcounts1)
dim(tddsNormcounts1)
rownames(tddsNormcounts1)
tddsNormcounts1 = data.frame(tddsNormcounts1)
write.table(tddsNormcounts1 , "tddsNormcounts1.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tddsNormcounts1)
tddsNormcounts1["Color"] <-  c("wt","wt","ko","ko")
dfx <-tddsNormcounts1[c(1:19995)]
PC<-prcomp(dfx)
head(PC)
PCi<-data.frame(PC$x,Color=tddsNormcounts1$Color)
percentage <- round(PC$sdev^2 / sum(PC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentage <- paste( colnames(PCi), "(", paste( as.character(percentage), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

p4<-ggplot(PCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentage[1]) + ylab(percentage[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("darkred","darkgreen","navy"))+
  scale_shape_manual(values=c(19,19,19,19,19,19))
p4 <- p4+theme_classic()
p4
ggsave("PCA_ddsNormcounts.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96)


#################################################### Transform data #####################################à
colnames(myTMM) <- c("wt_rep1", "wt_rep2", "ko_rep1", "ko_rep2")
head(myTMM)
boxplot(myTMM, ylim=c(0,500))
myTMM1 <- as.matrix(myTMM)
tmyTMM1 = t(myTMM1)
dim(tmyTMM1)
rownames(tmyTMM1)
tmyTMM1 = data.frame(tmyTMM1)
write.table(tmyTMM1 , "tmyTMM1.txt", sep="\t", quote = FALSE, append = FALSE)
dim(tmyTMM1)
tmyTMM1["Color"] <-  c("wt","wt","ko","ko")
dfx <-tmyTMM1[c(1:16581)]
PC<-prcomp(dfx)
head(PC)
PCi<-data.frame(PC$x,Color=tmyTMM1$Color)
percentage <- round(PC$sdev^2 / sum(PC$sdev^2) * 100, 2) #Plor variance var = sdev^2 https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
percentage <- paste( colnames(PCi), "(", paste( as.character(percentage), "%", ")", sep="") )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

p5<-ggplot(PCi,aes(x=PC1,y=PC2,col=Color, label=Color))+
  theme + xlab(percentage[1]) + ylab(percentage[2])+
  geom_point(size=2,alpha=1,aes(shape=Color))+
  scale_color_manual(values = c("darkred","darkgreen"))+
  scale_shape_manual(values=c(19,19,19,19,19,19))
p5 <- p5+theme_classic()
p5
ggsave("PCA_myTMM.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96)
#Volcano plot using ggplot2
p = ggplot(results, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(col=sig)) +
  scale_color_manual(values=c("red", "black"))+ xlim(-20,20)
p
#p+geom_text(data=filter(results, -log10(padj)>30), check_overlap = TRUE, size = 3, vjust = 0, hjust = 0, nudge_x = 0.05, aes(label=Genes))+theme_classic()
p+geom_text(data=results[(results$log2FoldChange > 5 | results$log2FoldChange < -5) & -log10(results$padj) >30,], check_overlap = TRUE, size = 3, vjust = 0, hjust = 0, nudge_x = 0.05 ,aes(label=Genes))+theme_classic()
icrgenes <- read.table("icr.genes",header=FALSE)
colnames(icrgenes) <- "id"
results_icr2 <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/results_icr2+.txt",header=TRUE)
head(results_icr2)
dim(results_icr2)
p = ggplot(results_icr2, aes(log2FoldChange, -log10(padj))) +
  geom_point()+
  geom_text(aes(label=ifelse(label != "no",as.character(Genes),'')),hjust=0,vjust=0)
p
p = ggplot(results_icr2, aes(log2FoldChange, -log10(padj), label =Genes)) +
  geom_point(aes(log2FoldChange, -log10(padj), color = results_icr2$label)) +  geom_label_repel(aes(label=ifelse(label != "no",as.character(Genes),'')),hjust=0,vjust=0,
                                                                                                box.padding   = 0.35, 
                                                                                                point.padding = 0.5,
                                                                                                segment.color = 'grey50') +  theme_classic()
p+scale_color_manual(values=c("red", "black"))
results_icr2 <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/results_icr2+.txt",header=TRUE)
head(results_icr2)
dim(results_icr2)
p = ggplot(results_icr2, aes(log2FoldChange, -log10(padj), label =Genes)) +
  geom_point(aes(log2FoldChange, -log10(padj), color = results_icr$label)) +  geom_label_repel(aes(label=ifelse(label != "no",as.character(Genes),'')),hjust=0,vjust=0,
                                                                                               box.padding   = 0.35, 
                                                                                               point.padding = 0.5,
                                                                                               segment.color = 'grey50') +  theme_classic()
p+scale_color_manual(values=c("red", "black"))

#EGSEA Pathway analysis
genes1 <- read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20_rearranged.txt", header = F)
head(genes1)
genes1 <- cSplit(genes1, "V4", ".")
head(genes1)
genes1 <- genes1$V4_1
DE <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/deseq2_results_res0.05_sorted.txt", header = T) 
head(DE)
DE <- data.frame(cbind(rownames(DE), DE))
head(DE)
colnames(DE) <- c("id",  "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "threshold")
logFC1 <- DE$log2FoldChange
DE <- cSplit(DE, "id", ".")
DE1 <- DE[,c(8,1:7)]
head(DE1)
dim(DE1)
DEgenes <- data.frame(DE1$id_1)
colnames(DEgenes) <- "id"
head(DEgenes)
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Mm.eg.db)
gene.df <- bitr(genes1, fromType = "ENSEMBL",
                toType = c( "ENTREZID", "SYMBOL"),
                OrgDb = org.Mm.eg.db)

dim(gene.df)
head(gene.df)
write.table(gene.df, "ensembletoentrez.id", row.names = F)
colnames(gene.df) <- c("id", "EntrezID", "Gene")
head(gene.df)
gene.df.id <- data.frame(gene.df$id)
colnames(gene.df.id) <- "id"
head(gene.df.id)
DEgenes_entrez = merge(DEgenes, gene.df, by="id", all.x=TRUE)
head(DEgenes_entrez)
dim(DEgenes_entrez)
write.table(DEgenes_entrez, "/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/DEgenes_entrez.txt", quote = F, append = F, row.names = F)

DEgenes_entrez_filt <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/DEgenes_entrez_filt.txt", header = F)
head(DEgenes_entrez_filt)
colnames(DEgenes_entrez_filt) <- c("id", "EntrezID", "logFC1")
dim(DEgenes_entrez_filt[,2:3])
DEgenes_entrez_symbol <- DEgenes_entrez_filt[,2:3]
logFC2 <- data.frame(cbind(DEgenes_entrez_filt$EntrezID, logFC1))
colnames(logFC2) <- c("EntrezID", "logFC1")
library(EGSEA)
gs.annots = buildIdx(entrezIDs = DEgenes_entrez_filt$EntrezID, species = "mouse",
                     gsdb.gsets = "all")
#gs.annots = buildIdx(entrezIDs = DEgenes_entrez_filt$EntrezID, species = "mouse",
#                     msigdb.gsets = "none", kegg.exclude = c("Metabolism"))
gsa = egsea.ora(geneIDs = DEgenes_entrez_filt$EntrezID, universe = as.character(gene.df$EntrezID), 
                logFC = logFC1, title = "WTvsZFP57KO", gs.annots = gs.annots,
                symbolsMap = DEgenes_entrez_symbol, display.top = 5, report.dir = "./WTvsZFP57KO",
                num.threads = 4, report = FALSE)
show(gsa)
names(gs.annots)
topSets(gsa, contrast = 1, gs.label = "kegg", number = 10)
t = topSets(gsa, contrast = 1, gs.label = "kegg", sort.by = "ora",
            number = 10, names.only = FALSE)
showSetByName(gsa, "c5", rownames(t)[1])
showSetByName(gsa, "kegg", rownames(t)[1])
plotMethods(gsa, gs.label = "kegg", contrast = 1, file.name = "WT-ZFP57KO-kegg-methods")
plotSummary(gsa, gs.label = "kegg", contrast = 1, file.name = "X24IL13-X24-kegg-summary")
showSetByID(gsa, gs.label = "kegg", c("hsa04060", "hsa04640"))
plotGOGraph(gsa, gs.label = "c5", file.name = "X24IL13-X24-c5-top-",
            sort.by = "avg.rank")
plotHeatmap(gsa, "Asthma", gs.label = "kegg", contrast = 1, file.name = "asthma-hm")
plotPathway(gsa, "Asthma", gs.label = "kegg", file.name = "asthma-pathway")
plotSummary(gsa, gs.label = "kegg", contrast = c(1, 2), file.name = "kegg-summary-cmp")
plotSummaryHeatmap(gsa, gs.label = "kegg", show.vals = "p.adj",
                   file.name = "il13-sum-heatmap")
plotHeatmap(gsa, "Asthma", gs.label = "kegg", contrast = "comparison",
            file.name = "asthma-hm-cmp")
plotPathway(gsa, "Asthma", gs.label = "kegg", contrast = 0, file.name = "asthma-pathway-cmp")


#Calculate FPKM
library("biomaRt")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
listFilters(mouse)
annot<-getBM(c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "strand", "start_position", "end_position","gene_biotype"), mart=mouse)
dim(annot)
annotsort <- annot[order(annot$ensembl_gene_id),]
head(annotsort)
annotsortLength <- data.frame(cbind(annotsort$ensembl_gene_id , annotsort$start_position, annotsort$end_position, annotsort$end_position - annotsort$start_position, annotsort$gene_biotype))
colnames(annotsortLength) <- c("Gene_ID","Start","End","length","Gene_type")
head(annotsortLength)
write.table(annotsortLength, "annotsortLength.txt", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE)
dim(annotsortLength)
annotsortLength.id <- data.frame(annotsortLength$Gene_ID)
write.table(annotsortLength.id, "annotsortLength.id", sep="\t", quote = FALSE, append = FALSE, row.names = FALSE, col.names = F)
head(annotsortLength.id) #56246
head(countdata)
dim(countdata)
mycounts_sorted <- countdata[order(rownames(countdata)),]
head(mycounts_sorted)
mycounts_sorted.id <- rownames(mycounts_sorted) #55604
write.table(mycounts_sorted, "mycounts_sorted.txt", sep="\t", quote = FALSE, append = FALSE)
write.table(mycounts_sorted.id, "mycounts_sorted.id", sep="\t", quote = FALSE, append = FALSE, row.names = F)
#Extract dotted length file (shared ids + 55274)
#fgrep -f annotsortLength.id mycounts_sorted.id | awk -F'.' '{print $1}' > annotsortLength_shared_mycounts_sorted.id
#fgrep -f annotsortLength.id mycounts_sorted.id -v | awk -F'.' '{print $1}' > annotsortLength_shared_mycounts_sortedNOT.id
#fgrep -f annotsortLength.id mycounts_sorted.id >  annotsortLength_shared_mycounts_sorted_dotted.id

#Extract overlapped  ID
#fgrep -f annotsortLength_shared_mycounts_sorted.id mycounts_sorted.txt > mycounts_sorted_refilt.txt
#fgrep -f annotsortLength_shared_mycounts_sorted.id annotsortLength.txt | sort -k1,1 -u > annotsortLength_refilt.txt
#paste annotsortLength_shared_mycounts_sorted_dotted.id annotsortLength_refilt.txt  | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6}' > annotsortLength_dotted.txt

mycounts_sorted_refilt <- read.table("mycounts_sorted_refilt.txt",header=FALSE)
rownames(mycounts_sorted_refilt)=mycounts_sorted_refilt[,1]
rownames(mycounts_sorted_refilt)
colnames(mycounts_sorted_refilt)
head(mycounts_sorted_refilt)
mycounts_sorted_refilt1 = mycounts_sorted_refilt[,-1]
mycounts_sorted_refilt1 = as.matrix(mycounts_sorted_refilt1)
head(mycounts_sorted_refilt1)
colnames(mycounts_sorted_refilt1) <- c("JB1_WT_Rep1","JB1_WT_Rep2","JB1_ZFP57_KO_Rep1","JB1_ZFP57_KO_Rep2")
head(mycounts_sorted_refilt1)
dim(mycounts_sorted_refilt1)

annotsortLength_dotted <- read.table("annotsortLength_dotted.txt",header=FALSE)
colnames(annotsortLength_dotted) <- c("Gene_ID","Start","End","length","Gene_type")
head(annotsortLength_dotted)
dim(annotsortLength_dotted)

samplematrix <- read.table("meanfragmentlength.txt",header=TRUE)
head(samplematrix)



featurelength <- read.table("featurecount_predicted_genelength_refilt.txt",header=FALSE)
featurelength <- featurelength$V2
head(featurelength)
dim(featurelength)
#DO NOT USE THIS WAY OF E-S GENE LENGTH:Gene End -start featurelength <- annotsortLength_dotted$length, only yhse featurecountpredicted
meanfragmentlength <- samplematrix$meanFragmentLength

#So now look commands
#counts = mycounts_sorted_refilt1
#gene.annotations = annotsortLength_dotted
#featureLength = featurelength
#samples.metrics = samplematrix
#meanFragmentLength = meanfragmentlength
# Return FPKM into a numeric matrix.
library(countToFPKM)
JB1_T12_fpkm_matrix <- fpkm (mycounts_sorted_refilt1, featurelength, meanfragmentlength)
head(JB1_T12_fpkm_matrix)
JB1_T12_fpkm_matrix <- data.frame(JB1_T12_fpkm_matrix)
JB1_T12_fpkm_matrix <- JB1_T12_fpkm_matrix[rowSums(JB1_T12_fpkm_matrix) > 0, ]

####write.table(JB1_T12_fpkm_matrix, "JB1_T12_fpkm_matrix_GeneLengthE_S.txt", sep="\t", quote = FALSE, append = FALSE, row.names = TRUE)
write.table(JB1_T12_fpkm_matrix, "JB1_T12_fpkm_matrix_GeneLengthFeatureCount.txt", sep="\t", quote = FALSE, append = FALSE, row.names = TRUE)

# Plot log10(FPKM+1) heatmap of top 30 highly variable features
fpkmheatmap(JB1_T12_fpkm_matrix, topvar=30, showfeaturenames=TRUE, return_log = TRUE)

#Calculate Log of FPKM
JB1_T12_fpkm_matrix_log <- log(JB1_T12_fpkm_matrix+1)
head(JB1_T12_fpkm_matrix_log)
write.table(JB1_T12_fpkm_matrix_log, "JB1_T12_fpkm_matrix_log_GeneLengthFeatureCount.txt", sep="\t", quote = FALSE, append = FALSE, row.names = TRUE)


#z-score calc on averaged normalized counts
#ddsNormcountsavg <- data.frame(cbind(rownames(ddsNormcounts),(ddsNormcounts$JB1_WT_Rep1+ddsNormcounts$JB1_WT_Rep2)/2, (ddsNormcounts$JB1_ZFP57_KO_Rep1 + ddsNormcounts$JB1_ZFP57_KO_Rep2)/2))
JB1_T12_fpkm_matrixavg <- data.frame(cbind((JB1_T12_fpkm_matrix$JB1_WT_Rep1+JB1_T12_fpkm_matrix$JB1_WT_Rep2)/2, (JB1_T12_fpkm_matrix$JB1_ZFP57_KO_Rep1 + JB1_T12_fpkm_matrix$JB1_ZFP57_KO_Rep2)/2))
rownames(JB1_T12_fpkm_matrixavg) <- rownames(JB1_T12_fpkm_matrix)
head(JB1_T12_fpkm_matrixavg)
dim(JB1_T12_fpkm_matrixavg)
colnames(JB1_T12_fpkm_matrixavg) <- c("JB1_WT","JB1_ZFP57KO")
rownames(JB1_T12_fpkm_matrixavg)
head(JB1_T12_fpkm_matrixavg)
####write.table(JB1_T12_fpkm_matrixavg, "JB1_T12_fpkm_matrix_GeneLengthE_Savg.txt", sep="\t", quote = FALSE, append = FALSE, row.names = TRUE)
write.table(JB1_T12_fpkm_matrixavg, "JB1_T12_fpkm_matrix_GeneLengthFeatureCountavg.txt", sep="\t", quote = FALSE, append = FALSE, row.names = TRUE)

#Plot by plot function
results_icr <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/results_icr.txt",header=TRUE)
head(results_icr)
dim(results_icr)
plot(results$log2FoldChange, -log10(results$padj), col="black", ylim=c(0,100), xlim=c(-20,20), pch=16, cex = 0.7, ylab="", xlab="", bty = 'n')
par(new=T)
plot(results_icr$log2FoldChange, -log10(results_icr$padj), col="blue", ylim=c(0,100), xlim=c(-20,20), pch=16 ,cex = 0.8, ylab="", xlab="", bty = 'n')

plot(results$log2FoldChange, -log10(results$padj), col="black", ylim=c(0,100), xlim=c(-10,10), pch=16, cex = 0.7, ylab="", xlab="", bty = 'n')
par(new=T)
plot(results_icr$log2FoldChange, -log10(results_icr$padj), col="blue", ylim=c(0,100), xlim=c(-10,10), pch=16 ,cex = 0.8, ylab="", xlab="", bty = 'n')
text(results_icr$log2FoldChange, -log10(results_icr$padj), results_icr$Genes , pos = 3, cex=0.6, col = "orange")


#Ratio based analysis
selected_markers_T0_T12_1 <- data.frame(selected_markers_T0_T12_1)
selected_markers_T0_T12_1.ratio <- cbind.data.frame(rownames(selected_markers_T0_T12_1),
                                                    (selected_markers_T0_T12_1$ZFP57KO_T0 / (selected_markers_T0_T12_1$WT_T0 + selected_markers_T0_T12_1$ZFP57KO_T0)), 
                                                    (selected_markers_T0_T12_1$ZFP57KO1_T12 / (selected_markers_T0_T12_1$WT1_T12 + selected_markers_T0_T12_1$ZFP57KO1_T12)),
                                                    (selected_markers_T0_T12_1$ZFP57KO2_T12 / (selected_markers_T0_T12_1$WT2_T12 + selected_markers_T0_T12_1$ZFP57KO2_T12)), stringsAsFactors=FALSE)
head(selected_markers_T0_T12_1.ratio)
rownames(selected_markers_T0_T12_1.ratio) <- selected_markers_T0_T12_1.ratio[,1]
selected_markers_T0_T12_1.ratio <- selected_markers_T0_T12_1.ratio[,-1]
head(selected_markers_T0_T12_1.ratio)
colnames(selected_markers_T0_T12_1.ratio) <- c("T0ratio", "T12ratio1", "T12ratio2")
head(selected_markers_T0_T12_1.ratio)
selected_markers_T0_T12_1.ratio1 <- data.frame(t(selected_markers_T0_T12_1.ratio))
head(selected_markers_T0_T12_1.ratio1)
selected_markers_T0_T12_1.ratio2 <- stack(selected_markers_T0_T12_1.ratio1)
head(selected_markers_T0_T12_1.ratio2)
selected_markers_T0_T12_1.ratio2["Group"] <- data.frame(rep(c("T0ratio", "T12ratio1", "T12ratio2"), 25)) #25 genes
head(selected_markers_T0_T12_1.ratio2)
selected_markers_T0_T12_1.ratio2 <- selected_markers_T0_T12_1.ratio2[,c(3,2,1)]
head(selected_markers_T0_T12_1.ratio2)
colnames(selected_markers_T0_T12_1.ratio2) <- c("Group", "Gene", "Ratio")
selected_markers_T0_T12_1.ratio2 <- data.frame(selected_markers_T0_T12_1.ratio2)
head(selected_markers_T0_T12_1.ratio2)
selected_markers_T0_T12_1.ratio2.0.5 <- selected_markers_T0_T12_1.ratio2
selected_markers_T0_T12_1.ratio2.0.5["Ratio"] <- selected_markers_T0_T12_1.ratio2$Ratio - 0.5

#Separate markers
p1 <-ggplot(data=selected_markers_T0_T12_1.ratio2.0.5[1:10,], aes(x=Gene, y=Ratio, label=Group, color=Group, fill=Group)) +
  geom_bar(width = 0.5,stat="identity",  position=position_dodge(), size=.3)+
  scale_fill_manual(values =(c("#009933", "#990099", "#990099"))) + scale_colour_manual(values =(c("#009933", "#990099", "#990099")))
p1 +  theme(panel.background = element_rect(fill = "#F5F5F5", size = 0.5, linetype = "solid"), axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(c(-1,1))

p2 <-ggplot(data=selected_markers_T0_T12_1.ratio2.0.5[11:24,], aes(x=Gene, y=Ratio, label=Group, color=Group, fill=Group)) +
  geom_bar(width = 0.5,stat="identity",  position=position_dodge(), size=.3)+
  scale_fill_manual(values =(c("#009933", "#990099", "#990099"))) + scale_colour_manual(values =(c("#009933", "#990099", "#990099")))
p2 +  theme(panel.background = element_rect(fill = "#F5F5F5", size = 0.5, linetype = "solid"), axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(c(-1,1))

p3 <-ggplot(data=selected_markers_T0_T12_1.ratio2.0.5[26:38,], aes(x=Gene, y=Ratio, label=Group, color=Group, fill=Group)) +
  geom_bar(width = 0.5,stat="identity",  position=position_dodge(), size=.3)+
  scale_fill_manual(values =(c("#009933", "#990099", "#990099"))) + scale_colour_manual(values =(c("#009933", "#990099", "#990099")))
p3 +  theme(panel.background = element_rect(fill = "#F5F5F5", size = 0.5, linetype = "solid"), axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(c(-1,1))

#Selected_markers_T0_T12_1Log <- log2(selected_markers_T0_T12_1+1)
colfunc <- colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","white","#FDDBC7","#F4A582","#D6604D","#B2182B"))
heatmap.2(selected_markers_T0_T12_1.ratio,trace = "none", col = colfunc , density.info=c("none"), dendrogram="none", scale = "none", sepwidth=c(0.001, 0.001), cexRow=1, font=3, cexCol = 1, margins =c(4,6), srtCol = 45, breaks = seq(0,10, length.out = 100))
#Dont sort order-already sorted
heatmap.2(selected_markers_T0_T12_1.ratio, 
          Colv = "NA", 
          Rowv = "NA", 
          trace = "none", 
          col = colfunc , 
          lmat=rbind( c(5, 4, 2), c(6,1,3)), 
          lwid=c(1, 4,2),
          lhei = c(1,7), 
          keysize=1.2, 
          key.par = list(cex=0.6), 
          density.info=c("none"), 
          dendrogram="none", 
          scale = "none", 
          sepwidth=c(0.001, 0.001), 
          cexRow=2, font=3, 
          cexCol = 1, 
          margins =c(5,8), 
          breaks = seq(0,1, length.out = 100))

# load package
library(pheatmap)
library(RColorBrewer)
breaksList = seq(0, 1, by = 0.1)
my_sample_row1 <- data.frame(Markers= c("ESC","ESC","ESC","ESC","ESC","Neural","Neural","Neural","Neural","Neural","Neural","Neural","Endo_Meso","Endo_Meso","Endo_Meso","Endo_Meso","Endo_Meso","Endo_Meso","Endo_Meso"))
rownames(my_sample_row1) <- rownames(selected_markers_T0_T12_1.ratio)
my_colour1 = list(Markers = c(ESC = "darkred", Neural = "orange", Endo_Meso = "magenta"))
#pheatmap(selected_markers_T0_T12_1.ratio,treeheight_row = 0, cluster_cols=F, cluster_rows=T, treeheight_col = 0, gaps_col =NULL, gaps_row = NULL, border_color = "black", breaks = breaksList,color= colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)))
pheatmap(selected_markers_T0_T12_1.ratio,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaksList)),
         breaks = breaksList,
         annotation_colors = my_colour1,
         annotation_row = my_sample_row1,
         cluster_cols=F, 
         cluster_rows=F,
         fontsize = 12,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "black",
         cutree_cols = 1)
#Calculate FPKM


#Calculate Log of avg FPKM
JB1_T12_fpkm_matrixavg_log <- log(data1+1)
head(JB1_T12_fpkm_matrixavg_log)
pheatmap(JB1_T12_fpkm_matrixavg_log,treeheight_row = 0, cluster_cols=F, cluster_rows=T, treeheight_col = 0, gaps_col =NULL, gaps_row = NULL, border_color = "black", breaks = breaksList,color= colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)))
write.table(JB1_T12_fpkm_matrixavg_log, "JB1_T12_fpkm_matrixavg_log_GeneLengthFeatureCount.txt", sep="\t", quote = FALSE, append = FALSE, row.names = TRUE)


library(gProfileR)
library(gprofiler2)
gene_annotation <- read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20_rearranged.txt", header = FALSE)
head(gene_annotation)
colnames(gene_annotation) <- c("chr","start","end","id","gene")
gene_annotation <- data.frame(gene_annotation)
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2")
deseq2_results_res0.05_gene_rearranged_filtered <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/all_deregulated_genes_deseq2_0.05.sorted_Zfp57_Kap1.txt", header = F)
head(deseq2_results_res0.05_gene_rearranged_filtered)
dim(deseq2_results_res0.05_gene_rearranged_filtered)
colnames(deseq2_results_res0.05_gene_rearranged_filtered) <- c("chr", "start", "end", "id", "Gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "threshold")
results_DE <- deseq2_results_res0.05_gene_rearranged_filtered
head(results_DE)
#Upregulated and Downregulated
results_DE.up <- results_DE[which(results_DE$log2FoldChange > 1),]
results_DE.up.id <- data.frame(results_DE.up$Gene)
colnames(results_DE.up.id) <- "id"
head(results_DE.up.id)
results_DE.up.symbol <- as.character(results_DE.up.id$id)

results_DE.down <- results_DE[which(results_DE$log2FoldChange < -1),]
results_DE.down.id <- data.frame(results_DE.down$Gene)
colnames(results_DE.down.id) <- "id"
results_DE.down.symbol <- as.character(results_DE.down.id$id)


results_DE.up.symbol <- read.table("/home/ankitv/Downloads/up100kb.txt")
results_DE.up.symbol <- data.frame(results_DE.up.symbol)
colnames(results_DE.down.id) <- "id"
head(results_DE.down.symbol)
results_DE.down.symbol <- as.character(results_DE.down.id)

gost.results_DE.up <- gost(query=results_DE.up.symbol, organism="mmusculus")
gostplot(gost.results_DE.up)
gost.results_DE.down <- gost(query=results_DE.down.symbol, organism="mmusculus")
gostplot(gost.results_DE.down)

#Sort by p-value
gost.results_DE.up.res.sorted <- gost.results_DE.up$result[order(gost.results_DE.up$result$p_value),]
head(gost.results_DE.up.res.sorted)
gost.results_DE.down.res.sorted <- gost.results_DE.down$result[order(gost.results_DE.down$result$p_value),]
head(gost.results_DE.down.res.sorted)

gost.results_DE.down.res.sorted_bar <- gost.results_DE.down.res.sorted[,c(11,3)]
head(gost.results_DE.down.res.sorted_bar)
gost.results_DE.down.res.sorted_bar_top <- head(gost.results_DE.down.res.sorted_bar,30)
gost.results_DE.down.res.sorted_bar_top$term_name <- gsub(' ', '.', gost.results_DE.down.res.sorted_bar_top$term_name)
gost.results_DE.down.res.sorted_bar_top$p_value <- -log10(gost.results_DE.down.res.sorted_bar_top$p_value)
gost.results_DE.down.res.sorted_bar_top <- data.frame(gost.results_DE.down.res.sorted_bar_top)

#rownames(gost.results_DE.down.res.sorted_bar_top) <- gost.results_DE.down.res.sorted_bar_top[,1]
#gost.results_DE.down.res.sorted_bar_top$term_name <- NULL
#par(las=2) # make label text perpendicular to axis
#par(mar=c(5,20,4,2)) # increase y-axis margin.
#barplot(-log10(gost.results_DE.down.res.sorted_bar_top$p_value), 
#        names.arg=rownames(gost.results_DE.down.res.sorted_bar_top), 
#        horiz=TRUE,  
#        cex.names=0.8,
#        col = "blue")
library(ggpubr)
ggbarplot(gost.results_DE.down.res.sorted_bar_top, x = "term_name", y = "p_value", color = "white", fill = "term_name" ,
          sort.by.groups = FALSE,
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "-log10(p.value)",
          xlab = "Pathways",
          legend.title = "Pathways",
          rotate = TRUE, position = position_dodge(),
          sort.val = "asc",
          ggtheme = theme_classic()
)
#gprofiler1
gPro.results_DE.up <- gprofiler(query=results_DE.up.symbol, organism="mmusculus", src_filter=NULL)
print(head(gPro.results_DE.up))
write.table(gPro.results_DE.up, "gPro.results_DE.up.txt", sep="\t", quote = FALSE, append = FALSE)
gPro.results_DE.down <- gprofiler(query=results_DE.down.symbol, organism="mmusculus", src_filter=NULL)
print(head(gPro.results_DE.down))
write.table(gPro.results_DE.down, "gPro.results_DE.down.txt", sep="\t", quote = FALSE, append = FALSE)


#Distance plot other way

#fgrep -f deseq2_results_res0.05_sorted.ensid /home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20_rearranged.txt -w | sort -k4,4 > deseq2_results_res0.05_sorted.ensid_chr.gene
#Merge using paste and add header and save as deseq2_results_res0.05_sorted.ensid_chr.gene_tagged.txt

#paste deseq2_results_res0.05_sorted.ensid_chr.gene deseq2_results_res0.05_sorted.txt > deseq2_results_res0.05_sorted.ensid_chr.gene_tagged.txt
#For interesection with bedtools: grep chr deseq2_results_res0.05_sorted.ensid_chr.gene_tagged.txt > all_deregulated_genes_deseq2_0.05.txt
#Filter non-windows/biased regions : For detailed description see Oct2019 in documents folder 
#sort -k1,1 -k2,2n all_deregulated_genes_deseq2_0.05.txt | grep ENS > all_deregulated_genes_deseq2_0.05.sorted1.txt
#bedtools intersect -wa -wb -a all_deregulated_genes_deseq2_0.05.sorted1.txt -b combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort_coverage_toremove.txt -v > all_deregulated_genes_deseq2_0.05.sorted.txt 
#164 imprinted_gene_name.txt
#fgrep -f imprinted_gene_name.txt all_deregulated_genes_deseq2_0.05.sorted.txt -w > imprinted_deregulated_genes_deseq2_0.05.sorted.txt 
#fgrep -f imprinted_gene_name.txt all_deregulated_genes_deseq2_0.05.sorted.txt -w -v > Nonimprinted_deregulated_genes_deseq2_0.05.sorted.txt

#Important do not remove duplicated genes because distance of genes may vary from one peak to another, sort low to high use sort -k9,9n
#All: bedtools closest -a all_deregulated_genes_deseq2_0.05.sorted.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8"\t"$12"\t"$17"\t"$19}' | sort -k9,9n > all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt
#Imprinted: bedtools closest -a imprinted_deregulated_genes_deseq2_0.05.sorted.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8"\t"$12"\t"$17"\t"$19}' | sort -k9,9n > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt
#Nonimprinted: bedtools closest -a Nonimprinted_deregulated_genes_deseq2_0.05.sorted.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8"\t"$12"\t"$17"\t"$19}' | sort -k9,9n > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt

#distanceplot <- read.table("distanceplot.txt",header=TRUE)
#head(distanceplot)
#plot(distanceplot$Group, distanceplot$All, type = "l", lwd = 3, main="Distance plot", col="black",xlim=c(0,9), ylim=c(0,100))
#par(new=T)
#plot(distanceplot$Group, distanceplot$Imprinted,  type = "l", lwd = 3, main="Distance plot", col="darkgreen",xlim=c(0,9), ylim=c(0,100))
#par(new=T)
#plot(distanceplot$Group, distanceplot$Nonimprinted,  type = "l", lwd = 3, main="Distance plot", col="#e56f12",xlim=c(0,9), ylim=c(0,100))


#p <-ggplot(data1, aes(x=WT, y=KO, color=Color)) + geom_point(aes(size =Color))+ 
scale_color_manual(values=c('red','grey', 'grey', 'grey','grey', 'grey','grey', 'grey'))+
  scale_size_manual(values=c(0.8,0.2,0.2,0.2,0.2,0.2,0.2,0.2))+
  geom_label_repel(aes(label=ifelse(Color == "DeregulatedImprintedUnder100Kb",as.character(Gene),'')),hjust=0,vjust=0,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')

#p + geom_abline(intercept = 0, color="black") +theme_classic()


#ggplot(data1, aes(x=WT, y=KO, color = Color)) + 
geom_point(aes(y = data1$Color, col = "y1")) + 
  geom_point(aes(y = data1$Gene, col = "y2"))

#Histogram of Zfp57 expression
hist(c(0.148085356,0.0535242681,0.0062007189,0.0038691993))

#Barplot of Zfp57 expression
rnazfp57bar <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/rnazfp57bar.txt",header = TRUE)
rownames(rnazfp57bar) <- rnazfp57bar$rep
rnazfp57bar["Group"] <- NULL
rnazfp57bar <- rnazfp57bar[,-1]
rnazfp57bar <- as.matrix(rnazfp57bar)
barplot(rnazfp57bar, col = c("green", "blue"))
rnazfp57bar1 <- stack(rnazfp57bar)
head(rnazfp57bar1)
colnames(rnazfp57bar1) <- c("exp", "sample")
library(ggpubr)
ggbarplot(rnazfp57bar1, x = "sample", y = "exp", color = "black", fill = "sample", sort.val = "asc", ylim = c(0,0.2),
          add = "mean_se", palette = c("darkgreen", "darkviolet"), legend.title = "Zfp57 RNA expression",
          position = position_dodge(), ggtheme = theme_classic())
ggsave("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/rnazfp57bar1.svg", width=12, height=8, units="cm", dpi=96)

#############################################################################################3

#Distance plot other way

#fgrep -f deseq2_results_res0.05_sorted.ensid /home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20_rearranged.txt -w | sort -k4,4 > deseq2_results_res0.05_sorted.ensid_chr.gene
#Merge using paste and add header and save as deseq2_results_res0.05_sorted.ensid_chr.gene_tagged.txt

#paste deseq2_results_res0.05_sorted.ensid_chr.gene deseq2_results_res0.05_sorted.txt > deseq2_results_res0.05_sorted.ensid_chr.gene_tagged.txt
#For interesection with bedtools: grep chr deseq2_results_res0.05_sorted.ensid_chr.gene_tagged.txt > all_deregulated_genes_deseq2_0.05.txt
#Filter non-windows/biased regions : For detailed description see Oct2019 in documents folder 
#sort -k1,1 -k2,2n all_deregulated_genes_deseq2_0.05.txt | grep ENS > all_deregulated_genes_deseq2_0.05.sorted1.txt
#bedtools intersect -wa -wb -a all_deregulated_genes_deseq2_0.05.sorted1.txt -b combined_NonwINDOWS_matRatioWTsort_matRatioZFP57KOsort_coverage_toremove.txt -v > all_deregulated_genes_deseq2_0.05.sorted.txt 
#164 imprinted_gene_name.txt
#fgrep -f imprinted_gene_name.txt all_deregulated_genes_deseq2_0.05.sorted.txt -w > imprinted_deregulated_genes_deseq2_0.05.sorted.txt 
#fgrep -f imprinted_gene_name.txt all_deregulated_genes_deseq2_0.05.sorted.txt -w -v > Nonimprinted_deregulated_genes_deseq2_0.05.sorted.txt

#Important do not remove duplicated genes because distance of genes may vary from one peak to another, sort low to high use sort -k9,9n
#All: bedtools closest -a all_deregulated_genes_deseq2_0.05.sorted.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8"\t"$12"\t"$17"\t"$19}' | sort -k9,9n > all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt
#Imprinted: bedtools closest -a imprinted_deregulated_genes_deseq2_0.05.sorted.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8"\t"$12"\t"$17"\t"$19}' | sort -k9,9n > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt
#Nonimprinted: bedtools closest -a Nonimprinted_deregulated_genes_deseq2_0.05.sorted.txt -b Zfp57_overlapped_Kap1_mm10_peaks.bed -d | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$8"\t"$12"\t"$17"\t"$19}' | sort -k9,9n > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt

#distanceplot <- read.table("distanceplot.txt",header=TRUE)
#head(distanceplot)
#plot(distanceplot$Group, distanceplot$All, type = "l", lwd = 3, main="Distance plot", col="black",xlim=c(0,9), ylim=c(0,100))
#par(new=T)
#plot(distanceplot$Group, distanceplot$Imprinted,  type = "l", lwd = 3, main="Distance plot", col="darkgreen",xlim=c(0,9), ylim=c(0,100))
#par(new=T)
#plot(distanceplot$Group, distanceplot$Nonimprinted,  type = "l", lwd = 3, main="Distance plot", col="#e56f12",xlim=c(0,9), ylim=c(0,100))


#p <-ggplot(data1, aes(x=WT, y=KO, color=Color)) + geom_point(aes(size =Color))+ 
scale_color_manual(values=c('red','grey', 'grey', 'grey','grey', 'grey','grey', 'grey'))+
  scale_size_manual(values=c(0.8,0.2,0.2,0.2,0.2,0.2,0.2,0.2))+
  geom_label_repel(aes(label=ifelse(Color == "DeregulatedImprintedUnder100Kb",as.character(Gene),'')),hjust=0,vjust=0,
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50')

#p + geom_abline(intercept = 0, color="black") +theme_classic()


#ggplot(data1, aes(x=WT, y=KO, color = Color)) + 
geom_point(aes(y = data1$Color, col = "y1")) + 
  geom_point(aes(y = data1$Gene, col = "y2"))

#Histogram of Zfp57 expression
hist(c(0.148085356,0.0535242681,0.0062007189,0.0038691993))

#Barplot of Zfp57 expression
rnazfp57bar <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/rnazfp57bar.txt",header = TRUE)
rownames(rnazfp57bar) <- rnazfp57bar$rep
rnazfp57bar["Group"] <- NULL
rnazfp57bar <- rnazfp57bar[,-1]
rnazfp57bar <- as.matrix(rnazfp57bar)
barplot(rnazfp57bar, col = c("green", "blue"))
rnazfp57bar1 <- stack(rnazfp57bar)
head(rnazfp57bar1)
colnames(rnazfp57bar1) <- c("exp", "sample")
library(ggpubr)
ggbarplot(rnazfp57bar1, x = "sample", y = "exp", color = "black", fill = "sample", sort.val = "asc", ylim = c(0,0.2),
          add = "mean_se", palette = c("darkgreen", "darkviolet"), legend.title = "Zfp57 RNA expression",
          position = position_dodge(), ggtheme = theme_classic())
ggsave("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/rnazfp57bar1.svg", width=12, height=8, units="cm", dpi=96)




#deseq2_mynoiseq_overlapped_0.05_padj.txt, something extra
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/")
data <- read.table("deseq2_mynoiseq_overlapped_0.05_padj.txt",header=TRUE)
rownames(data)
rownames(data)=data[,1]
rownames(data)
colnames(data)
head(data)
data1 = data.frame(cbind(data, -log10(data$padj)))
head(data1)
write.table(data1, "deseq2_mynoiseq_overlapped_0.05_padj_log.txt", sep="\t", quote = FALSE, append = FALSE, row.names = TRUE)








#Multiple Distance plot
fgrep -f set1gene.txt all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set1gene.txt
fgrep -f set1gene.txt imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set1gene.txt
fgrep -f set1gene.txt Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set1gene.txt
fgrep -f set2gene.txt all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set2gene.txt
fgrep -f set2gene.txt imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set2gene.txt
fgrep -f set2gene.txt Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set2gene.txt
fgrep -f set3gene.txt all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set3gene.txt
fgrep -f set3gene.txt imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set3gene.txt
fgrep -f set3gene.txt Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set3gene.txt
fgrep -f set4gene.txt all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set4gene.txt
fgrep -f set4gene.txt imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set4gene.txt
fgrep -f set4gene.txt Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set4gene.txt
fgrep -f set5gene.txt all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set5gene.txt
fgrep -f set5gene.txt imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set5gene.txt
fgrep -f set5gene.txt Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set5gene.txt
fgrep -f set6gene.txt all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set6gene.txt
fgrep -f set6gene.txt imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set6gene.txt
fgrep -f set6gene.txt Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set6gene.txt
fgrep -f set7gene.txt all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set7gene.txt
fgrep -f set7gene.txt imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set7gene.txt
fgrep -f set7gene.txt Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set7gene.txt
fgrep -f set8gene.txt all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set8gene.txt
fgrep -f set8gene.txt imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set8gene.txt
fgrep -f set8gene.txt Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set8gene.txt
fgrep -f set9gene.txt all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > all_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set9gene.txt
fgrep -f set9gene.txt imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set9gene.txt
fgrep -f set9gene.txt Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10.txt -w > Nonimprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_set9gene.txt

./distanceplot_awk_set1gene.sh > distanceplot_values_set1gene.txt
./distanceplot_awk_set2gene.sh > distanceplot_values_set2gene.txt
./distanceplot_awk_set3gene.sh > distanceplot_values_set3gene.txt
./distanceplot_awk_set4gene.sh > distanceplot_values_set4gene.txt
./distanceplot_awk_set5gene.sh > distanceplot_values_set5gene.txt
./distanceplot_awk_set6gene.sh > distanceplot_values_set6gene.txt
./distanceplot_awk_set7gene.sh > distanceplot_values_set7gene.txt
./distanceplot_awk_set8gene.sh > distanceplot_values_set8gene.txt
./distanceplot_awk_set9gene.sh > distanceplot_values_set9gene.txt




distanceplotvalueset1gene <- read.table("distanceplot_values_set1gene.txt", header=FALSE)
head(distanceplotvalueset1gene)
plotdistanceset1gene <- data.frame(c('0','10','100','1000','10000','100000','1000000','10000000','100000000','1000000000'))
plotdistanceset1gene["All"] <- as.numeric(as.character(distanceplotvalueset1gene[2:11,]))*100/as.numeric(as.character(distanceplotvalueset1gene[11,]))
plotdistanceset1gene["Imprinted"] <- as.numeric(as.character(distanceplotvalueset1gene[13:22,]))*100/as.numeric(as.character(distanceplotvalueset1gene[22,]))
plotdistanceset1gene["Nonmprinted"] <- as.numeric(as.character(distanceplotvalueset1gene[24:33,]))*100/as.numeric(as.character(distanceplotvalueset1gene[33,]))
colnames(plotdistanceset1gene) <- c("Group", "All", "Imprinted", "Nonimprinted")
plotdistanceset1gene
rownames(plotdistanceset1gene) <- plotdistanceset1gene$Group
plotdistanceset1gene["Group"] <- NULL
plotdistanceset1gene1 <- stack(plotdistanceset1gene)
head(plotdistanceset1gene1)
plotdistanceset1gene1["Distance"] <- c(0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9)
colnames(plotdistanceset1gene1) <- c("Genes", "Group", "Distance")
plotdistanceset1gene1 <- plotdistanceset1gene1[,c(2,3,1)]
head(plotdistanceset1gene1)
write.table(plotdistanceset1gene1,"/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/plotdistanceset1gene1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)
dset1gene <- ggplot(data=plotdistanceset1gene1, aes(x=Distance, y=Genes, group=Group)) +
  geom_line(aes(color=Group), size=1)+theme_classic()
dset1gene + scale_color_manual(values=c("black", "#00b050ff", "orange"))+ 
  scale_y_continuous(breaks = seq(0, 100, by = 10)) + 
  scale_x_continuous(breaks = seq(0, 100, by = 1))+
  labs(title="Distance plot", y ="% of deregulated genes in JB1", x = "Distance from JB1 ZFP57 peaks (bp)")+
  theme(axis.text.x = element_text(color="black", size=10, angle=315))
ggsave("distanceplotvalueset1gene.svg", width=12, height=8, units="cm", dpi=96)

distanceplotvalueset2gene <- read.table("distanceplot_values_set2gene.txt", header=FALSE)
head(distanceplotvalueset2gene)
plotdistanceset2gene <- data.frame(c('0','10','100','1000','10000','100000','1000000','10000000','100000000','1000000000'))
plotdistanceset2gene["All"] <- as.numeric(as.character(distanceplotvalueset2gene[2:11,]))*100/as.numeric(as.character(distanceplotvalueset2gene[11,]))
plotdistanceset2gene["Imprinted"] <- as.numeric(as.character(distanceplotvalueset2gene[13:22,]))*100/as.numeric(as.character(distanceplotvalueset2gene[22,]))
plotdistanceset2gene["Nonmprinted"] <- as.numeric(as.character(distanceplotvalueset2gene[24:33,]))*100/as.numeric(as.character(distanceplotvalueset2gene[33,]))
colnames(plotdistanceset2gene) <- c("Group", "All", "Imprinted", "Nonimprinted")
plotdistanceset2gene
rownames(plotdistanceset2gene) <- plotdistanceset2gene$Group
plotdistanceset2gene["Group"] <- NULL
plotdistanceset2gene1 <- stack(plotdistanceset2gene)
head(plotdistanceset2gene1)
plotdistanceset2gene1["Distance"] <- c(0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9)
colnames(plotdistanceset2gene1) <- c("Genes", "Group", "Distance")
plotdistanceset2gene1 <- plotdistanceset2gene1[,c(2,3,1)]
head(plotdistanceset2gene1)
write.table(plotdistanceset2gene1,"/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/plotdistanceset2gene1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)
dset2gene <- ggplot(data=plotdistanceset2gene1, aes(x=Distance, y=Genes, group=Group)) +
  geom_line(aes(color=Group), size=1)+theme_classic()
dset2gene + scale_color_manual(values=c("black", "#00b050ff", "orange"))+ 
  scale_y_continuous(breaks = seq(0, 100, by = 10)) + 
  scale_x_continuous(breaks = seq(0, 100, by = 1))+
  labs(title="Distance plot", y ="% of deregulated genes in JB1", x = "Distance from JB1 ZFP57 peaks (bp)")+
  theme(axis.text.x = element_text(color="black", size=10, angle=315))
ggsave("distanceplotvalueset2gene.svg", width=12, height=8, units="cm", dpi=96)


distanceplotvalueset3gene <- read.table("distanceplot_values_set3gene.txt", header=FALSE)
head(distanceplotvalueset3gene)
plotdistanceset3gene <- data.frame(c('0','10','100','1000','10000','100000','1000000','10000000','100000000','1000000000'))
plotdistanceset3gene["All"] <- as.numeric(as.character(distanceplotvalueset3gene[2:11,]))*100/as.numeric(as.character(distanceplotvalueset3gene[11,]))
plotdistanceset3gene["Imprinted"] <- as.numeric(as.character(distanceplotvalueset3gene[13:22,]))*100/as.numeric(as.character(distanceplotvalueset3gene[22,]))
plotdistanceset3gene["Nonmprinted"] <- as.numeric(as.character(distanceplotvalueset3gene[24:33,]))*100/as.numeric(as.character(distanceplotvalueset3gene[33,]))
colnames(plotdistanceset3gene) <- c("Group", "All", "Imprinted", "Nonimprinted")
plotdistanceset3gene
rownames(plotdistanceset3gene) <- plotdistanceset3gene$Group
plotdistanceset3gene["Group"] <- NULL
plotdistanceset3gene1 <- stack(plotdistanceset3gene)
head(plotdistanceset3gene1)
plotdistanceset3gene1["Distance"] <- c(0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9)
colnames(plotdistanceset3gene1) <- c("Genes", "Group", "Distance")
plotdistanceset3gene1 <- plotdistanceset3gene1[,c(2,3,1)]
head(plotdistanceset3gene1)
write.table(plotdistanceset3gene1,"/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/plotdistanceset3gene1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)
dset3gene <- ggplot(data=plotdistanceset3gene1, aes(x=Distance, y=Genes, group=Group)) +
  geom_line(aes(color=Group), size=1)+theme_classic()
dset3gene + scale_color_manual(values=c("black", "#00b050ff", "orange"))+ 
  scale_y_continuous(breaks = seq(0, 100, by = 10)) + 
  scale_x_continuous(breaks = seq(0, 100, by = 1))+
  labs(title="Distance plot", y ="% of deregulated genes in JB1", x = "Distance from JB1 ZFP57 peaks (bp)")+
  theme(axis.text.x = element_text(color="black", size=10, angle=315))
ggsave("distanceplotvalueset3gene.svg", width=12, height=8, units="cm", dpi=96)


distanceplotvalueset4gene <- read.table("distanceplot_values_set4gene.txt", header=FALSE)
head(distanceplotvalueset4gene)
plotdistanceset4gene <- data.frame(c('0','10','100','1000','10000','100000','1000000','10000000','100000000','1000000000'))
plotdistanceset4gene["All"] <- as.numeric(as.character(distanceplotvalueset4gene[2:11,]))*100/as.numeric(as.character(distanceplotvalueset4gene[11,]))
plotdistanceset4gene["Imprinted"] <- as.numeric(as.character(distanceplotvalueset4gene[13:22,]))*100/as.numeric(as.character(distanceplotvalueset4gene[22,]))
plotdistanceset4gene["Nonmprinted"] <- as.numeric(as.character(distanceplotvalueset4gene[24:33,]))*100/as.numeric(as.character(distanceplotvalueset4gene[33,]))
colnames(plotdistanceset4gene) <- c("Group", "All", "Imprinted", "Nonimprinted")
plotdistanceset4gene
rownames(plotdistanceset4gene) <- plotdistanceset4gene$Group
plotdistanceset4gene["Group"] <- NULL
plotdistanceset4gene1 <- stack(plotdistanceset4gene)
head(plotdistanceset4gene1)
plotdistanceset4gene1["Distance"] <- c(0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9)
colnames(plotdistanceset4gene1) <- c("Genes", "Group", "Distance")
plotdistanceset4gene1 <- plotdistanceset4gene1[,c(2,3,1)]
head(plotdistanceset4gene1)
write.table(plotdistanceset4gene1,"/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/plotdistanceset4gene1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)
dset4gene <- ggplot(data=plotdistanceset4gene1, aes(x=Distance, y=Genes, group=Group)) +
  geom_line(aes(color=Group), size=1)+theme_classic()
dset4gene + scale_color_manual(values=c("black", "#00b050ff", "orange"))+ 
  scale_y_continuous(breaks = seq(0, 100, by = 10)) + 
  scale_x_continuous(breaks = seq(0, 100, by = 1))+
  labs(title="Distance plot", y ="% of deregulated genes in JB1", x = "Distance from JB1 ZFP57 peaks (bp)")+
  theme(axis.text.x = element_text(color="black", size=10, angle=315))
ggsave("distanceplotvalueset4gene.svg", width=12, height=8, units="cm", dpi=96)


distanceplotvalueset5gene <- read.table("distanceplot_values_set5gene.txt", header=FALSE)
head(distanceplotvalueset5gene)
plotdistanceset5gene <- data.frame(c('0','10','100','1000','10000','100000','1000000','10000000','100000000','1000000000'))
plotdistanceset5gene["All"] <- as.numeric(as.character(distanceplotvalueset5gene[2:11,]))*100/as.numeric(as.character(distanceplotvalueset5gene[11,]))
plotdistanceset5gene["Imprinted"] <- as.numeric(as.character(distanceplotvalueset5gene[13:22,]))*100/as.numeric(as.character(distanceplotvalueset5gene[22,]))
plotdistanceset5gene["Nonmprinted"] <- as.numeric(as.character(distanceplotvalueset5gene[24:33,]))*100/as.numeric(as.character(distanceplotvalueset5gene[33,]))
colnames(plotdistanceset5gene) <- c("Group", "All", "Imprinted", "Nonimprinted")
plotdistanceset5gene
rownames(plotdistanceset5gene) <- plotdistanceset5gene$Group
plotdistanceset5gene["Group"] <- NULL
plotdistanceset5gene1 <- stack(plotdistanceset5gene)
head(plotdistanceset5gene1)
plotdistanceset5gene1["Distance"] <- c(0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9)
colnames(plotdistanceset5gene1) <- c("Genes", "Group", "Distance")
plotdistanceset5gene1 <- plotdistanceset5gene1[,c(2,3,1)]
head(plotdistanceset5gene1)
write.table(plotdistanceset5gene1,"/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/plotdistanceset5gene1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)
dset5gene <- ggplot(data=plotdistanceset5gene1, aes(x=Distance, y=Genes, group=Group)) +
  geom_line(aes(color=Group), size=1)+theme_classic()
dset5gene + scale_color_manual(values=c("black", "#00b050ff", "orange"))+ 
  scale_y_continuous(breaks = seq(0, 100, by = 10)) + 
  scale_x_continuous(breaks = seq(0, 100, by = 1))+
  labs(title="Distance plot", y ="% of deregulated genes in JB1", x = "Distance from JB1 ZFP57 peaks (bp)")+
  theme(axis.text.x = element_text(color="black", size=10, angle=315))
ggsave("distanceplotvalueset5gene.svg", width=12, height=8, units="cm", dpi=96)


distanceplotvalueset6gene <- read.table("distanceplot_values_set6gene.txt", header=FALSE)
head(distanceplotvalueset6gene)
plotdistanceset6gene <- data.frame(c('0','10','100','1000','10000','100000','1000000','10000000','100000000','1000000000'))
plotdistanceset6gene["All"] <- as.numeric(as.character(distanceplotvalueset6gene[2:11,]))*100/as.numeric(as.character(distanceplotvalueset6gene[11,]))
plotdistanceset6gene["Imprinted"] <- as.numeric(as.character(distanceplotvalueset6gene[13:22,]))*100/as.numeric(as.character(distanceplotvalueset6gene[22,]))
plotdistanceset6gene["Nonmprinted"] <- as.numeric(as.character(distanceplotvalueset6gene[24:33,]))*100/as.numeric(as.character(distanceplotvalueset6gene[33,]))
colnames(plotdistanceset6gene) <- c("Group", "All", "Imprinted", "Nonimprinted")
plotdistanceset6gene
rownames(plotdistanceset6gene) <- plotdistanceset6gene$Group
plotdistanceset6gene["Group"] <- NULL
plotdistanceset6gene1 <- stack(plotdistanceset6gene)
head(plotdistanceset6gene1)
plotdistanceset6gene1["Distance"] <- c(0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9)
colnames(plotdistanceset6gene1) <- c("Genes", "Group", "Distance")
plotdistanceset6gene1 <- plotdistanceset6gene1[,c(2,3,1)]
head(plotdistanceset6gene1)
write.table(plotdistanceset6gene1,"/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/plotdistanceset6gene1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)
dset6gene <- ggplot(data=plotdistanceset6gene1, aes(x=Distance, y=Genes, group=Group)) +
  geom_line(aes(color=Group), size=1)+theme_classic()
dset6gene + scale_color_manual(values=c("black", "#00b050ff", "orange"))+ 
  scale_y_continuous(breaks = seq(0, 100, by = 10)) + 
  scale_x_continuous(breaks = seq(0, 100, by = 1))+
  labs(title="Distance plot", y ="% of deregulated genes in JB1", x = "Distance from JB1 ZFP57 peaks (bp)")+
  theme(axis.text.x = element_text(color="black", size=10, angle=315))
ggsave("distanceplotvalueset6gene.svg", width=12, height=8, units="cm", dpi=96)


distanceplotvalueset7gene <- read.table("distanceplot_values_set7gene.txt", header=FALSE)
head(distanceplotvalueset7gene)
plotdistanceset7gene <- data.frame(c('0','10','100','1000','10000','100000','1000000','10000000','100000000','1000000000'))
plotdistanceset7gene["All"] <- as.numeric(as.character(distanceplotvalueset7gene[2:11,]))*100/as.numeric(as.character(distanceplotvalueset7gene[11,]))
plotdistanceset7gene["Imprinted"] <- as.numeric(as.character(distanceplotvalueset7gene[13:22,]))*100/as.numeric(as.character(distanceplotvalueset7gene[22,]))
plotdistanceset7gene["Nonmprinted"] <- as.numeric(as.character(distanceplotvalueset7gene[24:33,]))*100/as.numeric(as.character(distanceplotvalueset7gene[33,]))
colnames(plotdistanceset7gene) <- c("Group", "All", "Imprinted", "Nonimprinted")
plotdistanceset7gene
rownames(plotdistanceset7gene) <- plotdistanceset7gene$Group
plotdistanceset7gene["Group"] <- NULL
plotdistanceset7gene1 <- stack(plotdistanceset7gene)
head(plotdistanceset7gene1)
plotdistanceset7gene1["Distance"] <- c(0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9)
colnames(plotdistanceset7gene1) <- c("Genes", "Group", "Distance")
plotdistanceset7gene1 <- plotdistanceset7gene1[,c(2,3,1)]
head(plotdistanceset7gene1)
write.table(plotdistanceset7gene1,"/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/plotdistanceset7gene1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)
dset7gene <- ggplot(data=plotdistanceset7gene1, aes(x=Distance, y=Genes, group=Group)) +
  geom_line(aes(color=Group), size=1)+theme_classic()
dset7gene + scale_color_manual(values=c("black", "#00b050ff", "orange"))+ 
  scale_y_continuous(breaks = seq(0, 100, by = 10)) + 
  scale_x_continuous(breaks = seq(0, 100, by = 1))+
  labs(title="Distance plot", y ="% of deregulated genes in JB1", x = "Distance from JB1 ZFP57 peaks (bp)")+
  theme(axis.text.x = element_text(color="black", size=10, angle=315))
ggsave("distanceplotvalueset7gene.svg", width=12, height=8, units="cm", dpi=96)


distanceplotvalueset8gene <- read.table("distanceplot_values_set8gene.txt", header=FALSE)
head(distanceplotvalueset8gene)
plotdistanceset8gene <- data.frame(c('0','10','100','1000','10000','100000','1000000','10000000','100000000','1000000000'))
plotdistanceset8gene["All"] <- as.numeric(as.character(distanceplotvalueset8gene[2:11,]))*100/as.numeric(as.character(distanceplotvalueset8gene[11,]))
plotdistanceset8gene["Imprinted"] <- as.numeric(as.character(distanceplotvalueset8gene[13:22,]))*100/as.numeric(as.character(distanceplotvalueset8gene[22,]))
plotdistanceset8gene["Nonmprinted"] <- as.numeric(as.character(distanceplotvalueset8gene[24:33,]))*100/as.numeric(as.character(distanceplotvalueset8gene[33,]))
colnames(plotdistanceset8gene) <- c("Group", "All", "Imprinted", "Nonimprinted")
plotdistanceset8gene
rownames(plotdistanceset8gene) <- plotdistanceset8gene$Group
plotdistanceset8gene["Group"] <- NULL
plotdistanceset8gene1 <- stack(plotdistanceset8gene)
head(plotdistanceset8gene1)
plotdistanceset8gene1["Distance"] <- c(0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9)
colnames(plotdistanceset8gene1) <- c("Genes", "Group", "Distance")
plotdistanceset8gene1 <- plotdistanceset8gene1[,c(2,3,1)]
head(plotdistanceset8gene1)
write.table(plotdistanceset8gene1,"/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/plotdistanceset8gene1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)
dset8gene <- ggplot(data=plotdistanceset8gene1, aes(x=Distance, y=Genes, group=Group)) +
  geom_line(aes(color=Group), size=1)+theme_classic()
dset8gene + scale_color_manual(values=c("black", "#00b050ff", "orange"))+ 
  scale_y_continuous(breaks = seq(0, 100, by = 10)) + 
  scale_x_continuous(breaks = seq(0, 100, by = 1))+
  labs(title="Distance plot", y ="% of deregulated genes in JB1", x = "Distance from JB1 ZFP57 peaks (bp)")+
  theme(axis.text.x = element_text(color="black", size=10, angle=315))
ggsave("distanceplotvalueset8gene.svg", width=12, height=8, units="cm", dpi=96)


distanceplotvalueset9gene <- read.table("distanceplot_values_set9gene.txt", header=FALSE)
head(distanceplotvalueset9gene)
plotdistanceset9gene <- data.frame(c('0','10','100','1000','10000','100000','1000000','10000000','100000000','1000000000'))
plotdistanceset9gene["All"] <- as.numeric(as.character(distanceplotvalueset9gene[2:11,]))*100/as.numeric(as.character(distanceplotvalueset9gene[11,]))
plotdistanceset9gene["Imprinted"] <- as.numeric(as.character(distanceplotvalueset9gene[13:22,]))*100/as.numeric(as.character(distanceplotvalueset9gene[22,]))
plotdistanceset9gene["Nonmprinted"] <- as.numeric(as.character(distanceplotvalueset9gene[24:33,]))*100/as.numeric(as.character(distanceplotvalueset9gene[33,]))
colnames(plotdistanceset9gene) <- c("Group", "All", "Imprinted", "Nonimprinted")
plotdistanceset9gene
rownames(plotdistanceset9gene) <- plotdistanceset9gene$Group
plotdistanceset9gene["Group"] <- NULL
plotdistanceset9gene1 <- stack(plotdistanceset9gene)
head(plotdistanceset9gene1)
plotdistanceset9gene1["Distance"] <- c(0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8,9)
colnames(plotdistanceset9gene1) <- c("Genes", "Group", "Distance")
plotdistanceset9gene1 <- plotdistanceset9gene1[,c(2,3,1)]
head(plotdistanceset9gene1)
write.table(plotdistanceset9gene1,"/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/plotdistanceset9gene1.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)
dset9gene <- ggplot(data=plotdistanceset9gene1, aes(x=Distance, y=Genes, group=Group)) +
  geom_line(aes(color=Group), size=1)+theme_classic()
dset9gene + scale_color_manual(values=c("black", "#00b050ff", "orange"))+ 
  scale_y_continuous(breaks = seq(0, 100, by = 10)) + 
  scale_x_continuous(breaks = seq(0, 100, by = 1))+
  labs(title="Distance plot", y ="% of deregulated genes in JB1", x = "Distance from JB1 ZFP57 peaks (bp)")+
  theme(axis.text.x = element_text(color="black", size=10, angle=315))
ggsave("distanceplotvalueset9gene.svg", width=12, height=8, units="cm", dpi=96)

#Check alignment and all counts paste distanceplot_values* 


bedtools coverage -a ZFP57_bound_TE.txt -b ./../JB1_WT_Rep1_Aligned.sortedByCoord.out.bed > JB1_WT_Rep1_cov_Aligned.sortedByCoord.out.z.bed
bedtools coverage -a ZFP57_bound_TE.txt -b ./../JB1_WT_Rep2_Aligned.sortedByCoord.out.bed > JB1_WT_Rep2_cov_Aligned.sortedByCoord.out.z.bed
bedtools coverage -a ZFP57_bound_TE.txt -b ./../JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bed > JB1_ZFP57_KO_Rep1_cov_Aligned.sortedByCoord.out.z.bed
bedtools coverage -a ZFP57_bound_TE.txt -b ./../JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bed > JB1_ZFP57_KO_Rep2_cov_Aligned.sortedByCoord.out.z.bed

paste JB1_WT_Rep1_cov_Aligned.sortedByCoord.out.z.bed JB1_WT_Rep2_cov_Aligned.sortedByCoord.out.z.bed JB1_ZFP57_KO_Rep1_cov_Aligned.sortedByCoord.out.z.bed JB1_ZFP57_KO_Rep2_cov_Aligned.sortedByCoord.out.z.bed > mm10_ZFP57_bound_TE_bsplit.txt

setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2")
TEcovdata5 <- read.table("mm10_ZFP57_bound_TE_bsplit.txt", header=FALSE)
head(TEcovdata5)
tail(TEcovdata5)
dim(TEcovdata5)
#Sum B6 and JF1
TEcovdata5_reaaranged <- TEcovdata5[,c(1:5, 6, 15, 24, 33)]
head(TEcovdata5_reaaranged)
colnames(TEcovdata5_reaaranged) <- c("chr","start","end","strand","TE","JB1_WT_Rep1","JB1_WT_Rep2","JB1_ZFP57_KO_Rep1","JB1_ZFP57_KO_Rep2")
head(TEcovdata5_reaaranged)
rownames(TEcovdata5_reaaranged) <- paste0(TEcovdata5_reaaranged[,1], sep ="%", 
                                          TEcovdata5_reaaranged[,2], sep ="%",
                                          TEcovdata5_reaaranged[,3], sep ="%",
                                          TEcovdata5_reaaranged[,4], sep ="%",
                                          TEcovdata5_reaaranged[,5])
                                          
head(TEcovdata5_reaaranged)
colSums(TEcovdata5_reaaranged[,6:9])



TEcountData <- TEcovdata5_reaaranged[,6:9]


head(TEcountData)
TEcoldata <- read.table("coldata.txt" , header = TRUE, stringsAsFactors = FALSE)
rownames(TEcoldata) <- colnames(TEcountData)
head(TEcoldata)
TEcoldata <- TEcoldata[,c("condition","replicate")]
TEcoldata$condition <- factor(TEcoldata$condition)
TEcoldata$replicate <- factor(TEcoldata$replicate)
all(rownames(TEcoldata) == colnames(TEcountData)) #should print TRUE
TEdds <- DESeqDataSetFromMatrix(countData =TEcountData, colData = TEcoldata, design = ~ condition)
TEdds

keepTE <- rowSums(counts(TEdds)) >= 9
TEdds <- TEdds[keepTE,]
TEddscounts <- counts(TEdds, normalized=FALSE)
head(TEddscounts)
dim(TEddscounts)
colSums(TEddscounts)
#View filtered count matrix: View(counts(TEdds))
#Normalization is the part of DESeq command: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#Normalized separately:
TEddsNorm <- estimateSizeFactors(TEdds)
sizeFactors(TEddsNorm)
#Set size factor of bulk
sizeFactors(TEddsNorm) <- c(0.8706398,1.1225115,1.0022048,1.0605409)
sizeFactors(TEddsNorm)
#Note normalization =TRUE divide counts by the user set normalization factors
#Normalize allele specfic counts: export normalized counts
#Export normalized counts
TEddsNormcounts <- counts(TEddsNorm, normalized=TRUE)
head(TEddsNormcounts)
dim(TEddsNormcounts)
#write.table(TEddsNormcounts, "TEddsNormcounts.txt", sep="\t", quote=F, col.names=NA)
plotPCA(as.matrix(TEddsNormcounts)) #library(EDASeq)

#DEG
#Normalization is the part of DESeq command: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
TEdds <- DESeq(TEdds)
#See size factors estimated by DESeq function: sizeFactors(TEdds)
TEres <- results(TEdds, contrast=c("condition","ZFP57KO","WT"))
TEres
summary(TEres)
TEres_sorted = TEres[order(rownames(TEres)),]
head(TEres_sorted)
sum(TEres$padj < 0.05, na.rm=TRUE)
TEres_sorted = TEres[order(rownames(TEres)),]
head(TEres_sorted)
summary(TEres_sorted)
write.table(TEres_sorted, "deseq2_TEresults_TEres_sorted.txt", sep="\t", quote = FALSE, append = FALSE)
sum(TEres_sorted$padj < 0.05, na.rm=TRUE)
TEres_sorted$thTEreshold <- as.logical(TEres_sorted$padj < 0.05)
deseq2_TEresults_TEres0.05 <- TEres_sorted[which(TEres_sorted$thTEreshold == TRUE),]
head(deseq2_TEresults_TEres0.05)
dim(deseq2_TEresults_TEres0.05)
write.table(deseq2_TEresults_TEres0.05, "deseq2_TEresults_TEres0.05_sorted.txt", sep="\t", quote = FALSE, append = FALSE)

#NOISeq
myfactors = read.table("myfactors.txt", header=TRUE)
head(myfactors)
rownames(myfactors)
colnames(myfactors)

#Plot first PCA without filtering or normalisation of data: see ggplot2 
mydata1 <- NOISeq::readData(data=countforfpkm_filt1, length = gene.length, factors=myfactors)
myPCA1 = dat(mydata1, type = "PCA") #library(NOISeq)
par(mfrow = c(1, 2))
explo.plot(myPCA1, factor = "Sample")

#Normalisation
myRPKM = NOISeq::rpkm(assayData(mydata1)$exprs, long = gene.length, k = 0, lc = 1)
head(myRPKM)
mydata3 <- NOISeq::readData(data=myRPKM, factors=myfactors)
myPCA3 = dat(mydata3, type = "PCA")
par(mfrow = c(1, 2))
explo.plot(myPCA3, factor = "Sample")

#Low-count filtering
dim(myRPKM)
myRPKMFilt <- myRPKM[rowSums(myRPKM) > 0, ]
dim(myRPKMFilt)
myRPKMFilt_sorted = myRPKMFilt[order(rownames(myRPKMFilt)),]
head(myRPKMFilt_sorted)
dim(myRPKMFilt_sorted)
write.table(myRPKMFilt_sorted, "myRPKMFilt_sorted.txt", sep="\t", quote = FALSE, append = FALSE)
#Calculate Log of FPKM
myRPKMFilt_sorted_log <- log(myRPKMFilt_sorted+1)
head(myRPKMFilt_sorted_log)
write.table(myRPKMFilt_sorted_log, "myRPKMFilt_sorted_log_GeneLengthFeatureCount.txt", sep="\t", quote = FALSE, append = FALSE, row.names = TRUE)


#Motif analysis mm10 upper lower to upper
#Get fasta
bedtools getfasta -fi /home/ankitv/ref_av/mm10_bowtie2/mm10.fa -bed Zfp57_overlapped_Kap1_mm10_peaks.bed -fo Zfp57_overlapped_Kap1_mm10_peaks.fa

#Convert everything to upper case
awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' Zfp57_overlapped_Kap1_mm10_peaks.fa > Zfp57_overlapped_Kap1_mm10_peaks.upper.fa

#Run meme-chip
/home/ankitv/meme/bin/meme-chip -meme-nmotifs 2 Zfp57_overlapped_Kap1_mm10_peaks.upper.fa

bedtools getfasta -fi /home/ankitv/ref_av/mm10_bowtie2/mm10.fa -bed Zfp57_overlapped_Kap1_mm10_peaks.qualmore100.bed -fo Zfp57_overlapped_Kap1_mm10_peaks.qualmore100.fa

#Convert everything to upper case
awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' Zfp57_overlapped_Kap1_mm10_peaks.qualmore100.fa > Zfp57_overlapped_Kap1_mm10_peaks.qualmore100.upper.fa

#

bedtools getfasta -fi /home/ankitv/ref_av/mm10_bowtie2/mm10.fa -bed Zfp57_overlapped_Kap1_mm10_peaks.qualless100.bed -fo Zfp57_overlapped_Kap1_mm10_peaks.qualless100.fa

#Convert everything to upper case
awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' Zfp57_overlapped_Kap1_mm10_peaks.qualless100.fa > Zfp57_overlapped_Kap1_mm10_peaks.qualless100.upper.fa

#Assign Distance
bedtools closest -a Zfp57_overlapped_Kap1_mm10_peaks.bed -b /home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20_rearranged.txt -d > Zfp57_overlapped_Kap1_mm10_peaks_gene.bed


#Our Data Overlap with AF Smith Lab
#Peaks of AF Lab
cat GSE123942_ESC_ZFP57_PEAKS.bed  | sort -k1,1 -k2,2n > GSE123942_ALL_ZFP57_PEAKS.bed

bedtools intersect -wa -wb -a Zfp57_overlapped_Kap1_mm10_peaks.bed -b GSE123942_ALL_ZFP57_PEAKS.bed | awk '{print $4}' | sort -k1,1 -u  > Zfp57_overlapped_Kap1_mm10_GSE123942_ESC_merged.id
fgrep -f Zfp57_overlapped_Kap1_mm10_GSE123942_ESC_merged.id Zfp57_overlapped_Kap1_mm10_peaks.bed -w > Zfp57_overlapped_Kap1_mm10_peaks_overlapped_GSE123942_ESC.bed
bedtools closest -a Zfp57_overlapped_Kap1_mm10_peaks_overlapped_GSE123942_ESC.bed -b /home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20_rearranged.txt -d > Zfp57_overlapped_Kap1_mm10_peaks_overlapped_GSE123942_ESC.gene.bed

#Our Data Overlap with Anvar mm10 (Liftover.png of mm9 Zfp57peaks in excel sheet nar-02093-m-2015-File017.xlsx)
#Add peak column to Anvar_Zfp57_mm10.bed manually
bedtools intersect -wa -wb -a Zfp57_overlapped_Kap1_mm10_peaks.bed -b Anvar_Zfp57_mm10.bed | awk '{print $4}' | sort -k1,1 -u  > Zfp57_overlapped_Kap1_mm10_Anvar_Zfp57_mm10.id
fgrep -f Zfp57_overlapped_Kap1_mm10_Anvar_Zfp57_mm10.id Zfp57_overlapped_Kap1_mm10_peaks.bed -w > Zfp57_overlapped_Kap1_mm10_peaks_overlapped_Anvar_Zfp57.bed
bedtools closest -a Zfp57_overlapped_Kap1_mm10_peaks_overlapped_Anvar_Zfp57.bed -b /home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20_rearranged.txt -d > Zfp57_overlapped_Kap1_mm10_peaks_overlapped_Anvar_Zfp57.gene.bed

#Anvar mm10 Overlap with AF Smith Lab
bedtools intersect -wa -wb -a Anvar_Zfp57_mm10.bed -b GSE123942_ALL_ZFP57_PEAKS.bed | awk '{print $5}' | sort -k1,1 -u  > Anvar_Zfp57_mm10_GSE123942_ESC_merged.id
fgrep -f Anvar_Zfp57_mm10_GSE123942_ESC_merged.id Anvar_Zfp57_mm10.bed -w > Anvar_Zfp57_mm10_overlapped_GSE123942_ESC.bed

#Our Data Overlap with Anvar mm10 and then Overlap with AF Smith Lab
bedtools intersect -wa -wb -a Zfp57_overlapped_Kap1_mm10_peaks_overlapped_Anvar_Zfp57.bed -b GSE123942_ALL_ZFP57_PEAKS.bed | awk '{print $4}' | sort -k1,1 -u  > Zfp57_overlapped_Kap1_mm10_Anvar_Zfp57_mm10_GSE123942_ESC_merged.id


#GEO submition
#Give appropriate name to feature count matrix

#md5sum
c263d5a944619727eb39252784233e19  JB1_WT_Rep1_R1.fastq.gz
19409af90851a61e8cba00d693ee6b52  JB1_WT_Rep1_R2.fastq.gz
78ad5cf8dce4382d1dd7f73de343fd00  JB1_WT_Rep2_R1.fastq.gz
63114066beb94df00770d388e8832092  JB1_WT_Rep2_R2.fastq.gz
2f276109d3e3495048ed3d3612ca71da  JB1_ZFP57_KO_Rep1_R1.fastq.gz
02a3516d09684d88be25eb91f7789dfa  JB1_ZFP57_KO_Rep1_R2.fastq.gz
89cd6664924af095dfed294e09d3f9a7  JB1_ZFP57_KO_Rep2_R1.fastq.gz
f5fe953ee29d6b5ab800defa4414b1f5  JB1_ZFP57_KO_Rep2_R2.fastq.gz
