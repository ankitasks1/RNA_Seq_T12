#Bulk only hybrid
#folder in cluster 
#OurT0 hybrid JB1
/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotationlike-_chr.gtf -T 12 -o Bulk_hyb_ourT0_star-featureCounts_GRCm38.mm10.txt ourT0_WT_Rep1_Aligned.sortedByCoord.out.bam ourT0_WT_Rep2_Aligned.sortedByCoord.out.bam ourT0_WT_Rep3_Aligned.sortedByCoord.out.bam ourT0_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bam ourT0_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bam ourT0_ZFP57_KO_Rep3_Aligned.sortedByCoord.out.bam

#Our T0 Riso GSE77444, JB1
/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotationlike-_chr.gtf -T 12 -o Bulk_hyb_ourT0j_star-featureCounts_GRCm38.mm10.txt SRR3137618_Aligned.sortedByCoord.out.bam SRR3137619_Aligned.sortedByCoord.out.bam

#Bouschet
/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotationlike-_chr.gtf -T 12 -o Bulk_hyb_bouchet_star-featureCounts_GRCm38.mm10.txt SRR1409932_Aligned.sortedByCoord.out.bam SRR1409933_Aligned.sortedByCoord.out.bam SRR1409934_Aligned.sortedByCoord.out.bam SRR1409935_Aligned.sortedByCoord.out.bam SRR1409938_Aligned.sortedByCoord.out.bam SRR1409939_Aligned.sortedByCoord.out.bam SRR1409940_Aligned.sortedByCoord.out.bam SRR1409941_Aligned.sortedByCoord.out.bam

#T12
/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -p -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotationlike-_chr.gtf -T 12 -o Bulk_hyb_ourT12_star-featureCounts_GRCm38.mm10.txt ourT12_WT_Rep1_Aligned.sortedByCoord.out.bam ourT12_WT_Rep2_Aligned.sortedByCoord.out.bam ourT12_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bam ourT12_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bam

scp ankitv@lilligridbio.na.iac.cnr.it:/home/ankitv/rna-seq/mouse/2020/combine_ourn_bouschet/Bulk_hyb_ourT0_star-featureCounts_GRCm38.mm10.txt ./
scp ankitv@lilligridbio.na.iac.cnr.it:/home/ankitv/rna-seq/mouse/2020/combine_ourn_bouschet/Bulk_hyb_ourT0j_star-featureCounts_GRCm38.mm10.txt ./
scp ankitv@lilligridbio.na.iac.cnr.it:/home/ankitv/rna-seq/mouse/2020/combine_ourn_bouschet/Bulk_hyb_bouchet_star-featureCounts_GRCm38.mm10.txt ./
scp ankitv@lilligridbio.na.iac.cnr.it:/home/ankitv/rna-seq/mouse/2020/combine_ourn_bouschet/Bulk_hyb_ourT12_star-featureCounts_GRCm38.mm10.txt ./
  
  
#Bulk only hybrid
  
################################################################# DESeq2 ######################################################################
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2020/CombineOurNBouschetNFeilNFerguson/")
bhybcountdata1 <- read.table("Bulk_hyb_ourT0j_star-featureCounts_GRCm38.mm10.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
bhybcountdata1 <- bhybcountdata1[ ,6:ncol(bhybcountdata1)]
head(bhybcountdata1)
# Remove .bam or .sam from filenames
colnames(bhybcountdata1) <- gsub("\\.[sb]am$", "", colnames(bhybcountdata1))
colnames(bhybcountdata1) <-  c("ourT0WTJB1","ourT0ZFP57KOJB1")
bhybcountdata1 <- bhybcountdata1[order(rownames(bhybcountdata1)),]
head(bhybcountdata1)
dim(bhybcountdata1)
tail(bhybcountdata1)
library(NOISeq)
bhybcountdata1uq <- uqua(bhybcountdata1, long = 1000, lc = 0, k = 0)
head(bhybcountdata1uq)


bhybcountdata2 <- read.table("Bulk_hyb_bouchet_star-featureCounts_GRCm38.mm10.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
bhybcountdata2 <- bhybcountdata2[ ,6:ncol(bhybcountdata2)]
head(bhybcountdata2)
# Remove .bam or .sam from filenames
colnames(bhybcountdata2) <- gsub("\\.[sb]am$", "", colnames(bhybcountdata2))
colnames(bhybcountdata2) <-  c("SRR1409932","SRR1409933","SRR1409934","SRR1409935","SRR1409938","SRR1409939","SRR1409940","SRR1409941")
bhybcountdata2 <- bhybcountdata2[order(rownames(bhybcountdata2)),]
head(bhybcountdata2)
dim(bhybcountdata2)
tail(bhybcountdata2)
bhybcountdata2uq <- uqua(bhybcountdata2, long = 1000, lc = 0, k = 0)
head(bhybcountdata2uq)

bhybcountdata3 <- read.table("Bulk_hyb_ourT12_star-featureCounts_GRCm38.mm10.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
bhybcountdata3 <- bhybcountdata3[ ,6:ncol(bhybcountdata3)]
head(bhybcountdata3)
# Remove .bam or .sam from filenames
colnames(bhybcountdata3) <- gsub("\\.[sb]am$", "", colnames(bhybcountdata3))
colnames(bhybcountdata3) <-  c("ourT12WT1","ourT12WT2","ourT12ZFP57KO1","ourT12ZFP57KO2")
bhybcountdata3 <- bhybcountdata3[order(rownames(bhybcountdata3)),]
head(bhybcountdata3)
dim(bhybcountdata3)
tail(bhybcountdata3)
bhybcountdata3uq <- uqua(bhybcountdata3, long = 1000, lc = 0, k = 0)
head(bhybcountdata3uq)


bhybcountdata <- cbind.data.frame(rownames(bhybcountdata1),bhybcountdata1, rownames(bhybcountdata2), bhybcountdata2, rownames(bhybcountdata3), bhybcountdata3)
head(bhybcountdata)
tail(bhybcountdata)
dim(bhybcountdata)
bhybcountdata <- bhybcountdata[,c(2,3,5:12,14:17)]
head(bhybcountdata)
library(DESeq2)
bhybcountdata = as.matrix(bhybcountdata)
head(bhybcountdata)
dim(bhybcountdata)
boxplot(bhybcountdata, ylim=c(0,100))



bhybcountdatauq <- cbind.data.frame(bhybcountdata1uq, bhybcountdata2uq, bhybcountdata3uq)
bhybcountdatauq <- cbind.data.frame(rownames(bhybcountdata1uq),bhybcountdata1uq, rownames(bhybcountdata2uq), bhybcountdata2uq, rownames(bhybcountdata3uq), bhybcountdata3uq)
head(bhybcountdatauq)
tail(bhybcountdatauq)
dim(bhybcountdatauq)
bhybcountdatauq <- bhybcountdatauq[,c(2,3,5:12,14:17)]
library(DESeq2)
bhybcountdatauq = data.frame(bhybcountdatauq)
head(bhybcountdatauq)
dim(bhybcountdatauq)
plotPCA(as.matrix(bhybcountdatauq[,1:14]))
#boxplot(bhybcountdatauq, ylim=c(0,100))
bhybcountdatauq["id"] <- rownames(bhybcountdatauq)


#created after basila gene selection  bouschetmarkerslist.sorted.txt + /media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/selected.markers
markers.gene.list.new <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/markers.gene.list.new.txt", header = F)
head(markers.gene.list.new)
colnames(markers.gene.list.new) <- c("Gene", "markers")
dim(markers.gene.list.new) #print 41 genes
ens_gene_names_chrpos_dedup_M20 <- read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt")
head(ens_gene_names_chrpos_dedup_M20)
dim(ens_gene_names_chrpos_dedup_M20)
ens_gene_names_chrpos_dedup_M20 <- ens_gene_names_chrpos_dedup_M20[,c(3,4,5,1,2)]
head(ens_gene_names_chrpos_dedup_M20)
colnames(ens_gene_names_chrpos_dedup_M20) <- c("chr", "start", "end", "id", "Gene")

cellmarkers = merge(markers.gene.list.new, ens_gene_names_chrpos_dedup_M20, by="Gene", all.x=FALSE, sort = FALSE)
head(cellmarkers)
dim(cellmarkers) #print all 41 genes assigned


bhybcountdatauq.mergemarkers = merge(cellmarkers,bhybcountdatauq, by="id", all.x=FALSE, sort = FALSE)
head(bhybcountdatauq.mergemarkers)
dim(bhybcountdatauq.mergemarkers) #Print all 41 genes overlapped and assigned
bhybcountdatauq.cellmarkers <- bhybcountdatauq.mergemarkers
dim(bhybcountdatauq.cellmarkers)
head(bhybcountdatauq.cellmarkers)
bhybcountdatauq.cellmarkersre <- bhybcountdatauq.cellmarkers[,c(2, 7:20)]
head(bhybcountdatauq.cellmarkersre,2)
rownames(bhybcountdatauq.cellmarkersre) <- bhybcountdatauq.cellmarkersre[,1]
head(bhybcountdatauq.cellmarkersre)
bhybcountdatauq.cellmarkers1 <- bhybcountdatauq.cellmarkersre[,c(2:15)]
head(bhybcountdatauq.cellmarkers1,2)
dim(bhybcountdatauq.cellmarkers1)
bhybcountdatauq.cellmarkers1 <- as.matrix(bhybcountdatauq.cellmarkers1)

#Loading control
loadingcontrol.cords <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/loadingcontrol.symbol.ens.chr.txt", header = F)
head(loadingcontrol.cords)
colnames(loadingcontrol.cords) <- c("chr", "start", "end", "id", "Gene")
bhybcountdatauq.loadingcontrol = merge(loadingcontrol.cords,bhybcountdatauq, by="id", all.x=FALSE)
head(bhybcountdatauq.loadingcontrol)
dim(bhybcountdatauq.loadingcontrol)
bhybcountdatauq.loadingcontrol <- bhybcountdatauq.loadingcontrol[,c(5:19)]
head(bhybcountdatauq.loadingcontrol,2)
rownames(bhybcountdatauq.loadingcontrol) <- bhybcountdatauq.loadingcontrol$Gene
bhybcountdatauq.loadingcontrol1 <- bhybcountdatauq.loadingcontrol[,c(2:15)]
head(bhybcountdatauq.loadingcontrol1,2)
bhybcountdatauq.loadingcontrol1 <- as.matrix(bhybcountdatauq.loadingcontrol1)
bhybcountdatauq.loadingcontrol1.Actb <- bhybcountdatauq.loadingcontrol1[1,]
#Normalize with loading control
dim(bhybcountdatauq.cellmarkers1)

bhybcountdatauq.cellmarkers1.norm <- cbind.data.frame((bhybcountdatauq.cellmarkers1[,1] / bhybcountdatauq.loadingcontrol1.Actb[1]),
                                                      (bhybcountdatauq.cellmarkers1[,2] / bhybcountdatauq.loadingcontrol1.Actb[2]),
                                                      (bhybcountdatauq.cellmarkers1[,3] / bhybcountdatauq.loadingcontrol1.Actb[3]),
                                                      (bhybcountdatauq.cellmarkers1[,4] / bhybcountdatauq.loadingcontrol1.Actb[4]),
                                                      (bhybcountdatauq.cellmarkers1[,5] / bhybcountdatauq.loadingcontrol1.Actb[5]),
                                                      (bhybcountdatauq.cellmarkers1[,6] / bhybcountdatauq.loadingcontrol1.Actb[6]),
                                                      (bhybcountdatauq.cellmarkers1[,7] / bhybcountdatauq.loadingcontrol1.Actb[7]),
                                                      (bhybcountdatauq.cellmarkers1[,8] / bhybcountdatauq.loadingcontrol1.Actb[8]),
                                                      (bhybcountdatauq.cellmarkers1[,9] / bhybcountdatauq.loadingcontrol1.Actb[9]),
                                                      (bhybcountdatauq.cellmarkers1[,10] / bhybcountdatauq.loadingcontrol1.Actb[10]),
                                                      (bhybcountdatauq.cellmarkers1[,11] / bhybcountdatauq.loadingcontrol1.Actb[11]),
                                                      (bhybcountdatauq.cellmarkers1[,12] / bhybcountdatauq.loadingcontrol1.Actb[12]),
                                                      (bhybcountdatauq.cellmarkers1[,13] / bhybcountdatauq.loadingcontrol1.Actb[13]),
                                                      (bhybcountdatauq.cellmarkers1[,14] / bhybcountdatauq.loadingcontrol1.Actb[14]))
head(bhybcountdatauq.cellmarkers1.norm)
colnames(bhybcountdatauq.cellmarkers1.norm) <- c("ourT0WTJB1","ourT0ZFP57KOJB1","SRR1409932","SRR1409933","SRR1409934","SRR1409935","SRR1409938","SRR1409939","SRR1409940","SRR1409941","ourT12WT1", "ourT12WT2", "ourT12ZFP57KO1", "ourT12ZFP57KO2")

head(bhybcountdatauq.cellmarkers1.norm)
boxplot(log10(bhybcountdatauq.cellmarkers1.norm), ylim=c(0,10), col = c("blue","blue","green","green","green","green","green","green","green","green","orange","orange","orange","orange"))
#save as bhybcountdatauq.cellmarkers1.norm.sep.svg
library(EDASeq)
plotPCA(as.matrix(bhybcountdatauq.cellmarkers1.norm))
plotPCA(as.matrix(bhybcountdatauq.cellmarkers1.norm), labels=F, col =  c("cyan","red","blue","blue","darkgreen","darkgreen","blue","blue","darkgreen","darkgreen","navy","navy","red","red"))
ggsave("PCA_tbhybcountdatauq.cellmarkers1.norm_scaleT.sep.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

write.table(bhybcountdatauq.cellmarkers1.norm, "bhybcountdatauq.cellmarkers1.norm.sep.txt", row.names = T, quote = F, append=F)
#Dendrogram cluster
dim(bhybcountdatauq.cellmarkers1.norm)  #41 14
# Compute distances and hierarchical clustering
bhybcountdatauq.cellmarkers1.norm1 <- t(bhybcountdatauq.cellmarkers1.norm)
bhybcountdatauq.cellmarkers1.normscale <- scale(bhybcountdatauq.cellmarkers1.norm1)
ddbhybcountdatauq.cellmarkers1.normscale <- dist(as.matrix(bhybcountdatauq.cellmarkers1.normscale), method = "euclidean")
ddbhybcountdatauq.cellmarkers1.normscalehc <- hclust(ddbhybcountdatauq.cellmarkers1.normscale, method = "complete")
plot(ddbhybcountdatauq.cellmarkers1.normscalehc)
#save as ddbhybcountdatauq.cellmarkers1.normscalehc.sep.svg
head(bhybcountdatauq.cellmarkers1.norm,1)                                                             
colnames(bhybcountdatauq.cellmarkers1.norm) <- colnames(bhybcountdatauq.cellmarkers1)
head(bhybcountdatauq.cellmarkers1.norm, 2)
colfunc <- colorRampPalette(c("navy","white","firebrick3"))
head(bhybcountdatauq.cellmarkers1.norm)
dim(bhybcountdatauq.cellmarkers1.norm)  #41 14
summary(bhybcountdatauq.cellmarkers1.norm)
library(pheatmap)
library(RColorBrewer)
breaksList2 = seq(0, 0.05, by = 0.001)
my_sample_row1 <- data.frame(Markers= c("ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Knockout"))
rownames(my_sample_row1) <- rownames(bhybcountdatauq.cellmarkers1.norm)
my_colour1 = list(Markers = c(ESC = "#00ff00", Neural.cortex ="#008000", Endo.Meso = "#ff80b3", Knockout = "#B3A580"))
#pheatmap(selected_markers_T0_T12_1.ratio,treeheight_row = 0, cluster_cols=F, cluster_rows=T, treeheight_col = 0, gaps_col =NULL, gaps_row = NULL, border_color = "black", breaks = breaksList,color= colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)))
#Remove Otx2
pheatmap(bhybcountdatauq.cellmarkers1.norm[c(1:17,19:41),],
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(length(breaksList2)),
         breaks = breaksList2,
         treeheight_row = 30,
         treeheight_col = 30,
         annotation_colors = my_colour1,
         fontsize = 8,
         annotation_row = my_sample_row1,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=T, 
         cluster_rows=F,
         cutree_cols = 3,
         cellwidth = 10, 
         cellheight = 10)
#save as filename = "pheatmap_bhybcountdatauq.cellmarkers1.norm.sep.png")

dev.off()

#Z-score
z_Tbhybcountdatauq.cellmarkers1.norm= scale(t(bhybcountdatauq.cellmarkers1.norm), center = TRUE, scale = TRUE)
z_bhybcountdatauq.cellmarkers1.norm <- t(z_Tbhybcountdatauq.cellmarkers1.norm)
colfunc <- colorRampPalette(c("navy","white","firebrick3"))
head(z_bhybcountdatauq.cellmarkers1.norm)
dim(z_bhybcountdatauq.cellmarkers1.norm)  #41 14
summary(z_bhybcountdatauq.cellmarkers1.norm)
library(pheatmap)
library(RColorBrewer)
breaksList1 = seq(-1, 1, by = 0.01)
my_sample_row1 <- data.frame(Markers= c("ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Knockout"))
rownames(my_sample_row1) <- rownames(z_bhybcountdatauq.cellmarkers1.norm)
my_colour1 = list(Markers = c(ESC = "#ebeb14", Neural.cortex ="#8ca50e", Endo.Meso = "#ff80b3", Knockout = "#B3A580"))
#pheatmap(selected_markers_T0_T12_1.ratio,treeheight_row = 0, cluster_cols=F, cluster_rows=T, treeheight_col = 0, gaps_col =NULL, gaps_row = NULL, border_color = "black", breaks = breaksList,color= colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)))
#Remove Otx2
pheatmap(z_bhybcountdatauq.cellmarkers1.norm[c(1:17,19:41),],
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(length(breaksList1)),
         breaks = breaksList1,
         treeheight_row = 30,
         treeheight_col = 30,
         annotation_colors = my_colour1,
         fontsize = 8,
         annotation_row = my_sample_row1,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=T, 
         cluster_rows=F,
         cutree_cols = 2,
         cellwidth = 10,
         cellheight = 10)

#save as filename = "pheatmap_z_bhybcountdatauq.cellmarkers1.norm.sep.svg")

dev.off()

#xx
bhybcountdatauqre <- bhybcountdatauq
head(bhybcountdatauqre,2)
dim(bhybcountdatauqre)
#Adjust by Actb
bhybcountdatauqre.norm <- cbind.data.frame((rownames(bhybcountdatauqre)),
                                           (bhybcountdatauqre[,1] / bhybcountdatauq.loadingcontrol1.Actb[1]),
                                           (bhybcountdatauqre[,2] / bhybcountdatauq.loadingcontrol1.Actb[2]),
                                           (bhybcountdatauqre[,3] / bhybcountdatauq.loadingcontrol1.Actb[3]),
                                           (bhybcountdatauqre[,4] / bhybcountdatauq.loadingcontrol1.Actb[4]),
                                           (bhybcountdatauqre[,5] / bhybcountdatauq.loadingcontrol1.Actb[5]),
                                           (bhybcountdatauqre[,6] / bhybcountdatauq.loadingcontrol1.Actb[6]),
                                           (bhybcountdatauqre[,7] / bhybcountdatauq.loadingcontrol1.Actb[7]),
                                           (bhybcountdatauqre[,8] / bhybcountdatauq.loadingcontrol1.Actb[8]),
                                           (bhybcountdatauqre[,9] / bhybcountdatauq.loadingcontrol1.Actb[9]),
                                           (bhybcountdatauqre[,10] / bhybcountdatauq.loadingcontrol1.Actb[10]),
                                           (bhybcountdatauqre[,11] / bhybcountdatauq.loadingcontrol1.Actb[11]),
                                           (bhybcountdatauqre[,12] / bhybcountdatauq.loadingcontrol1.Actb[12]),
                                           (bhybcountdatauqre[,13] / bhybcountdatauq.loadingcontrol1.Actb[13]),
                                           (bhybcountdatauqre[,14] / bhybcountdatauq.loadingcontrol1.Actb[14]))
head(bhybcountdatauqre.norm)
colnames(bhybcountdatauqre.norm) <- c("id","ourT0WTJB1","ourT0ZFP57KOJB1","SRR1409932","SRR1409933","SRR1409934","SRR1409935","SRR1409938","SRR1409939","SRR1409940","SRR1409941","ourT12WT1", "ourT12WT2", "ourT12ZFP57KO1", "ourT12ZFP57KO2")

head(bhybcountdatauqre.norm)
rownames(bhybcountdatauqre.norm) <- bhybcountdatauqre.norm$id

head(bhybcountdatauqre.norm,2)
dim(bhybcountdatauqre.norm)

#Extract cell marker coordinates
#created after basila gene selection  bouschetmarkerslist.sorted.txt + /media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/selected.markers
markers.gene.list.new <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/markers.gene.list.new.txt", header = F)
head(markers.gene.list.new)
colnames(markers.gene.list.new) <- c("Gene", "markers")
dim(markers.gene.list.new) #print 41 genes
ens_gene_names_chrpos_dedup_M20 <- read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt")
head(ens_gene_names_chrpos_dedup_M20)
dim(ens_gene_names_chrpos_dedup_M20)
ens_gene_names_chrpos_dedup_M20 <- ens_gene_names_chrpos_dedup_M20[,c(3,4,5,1,2)]
head(ens_gene_names_chrpos_dedup_M20)
colnames(ens_gene_names_chrpos_dedup_M20) <- c("chr", "start", "end", "id", "Gene")

cellmarkers = merge(markers.gene.list.new, ens_gene_names_chrpos_dedup_M20, by="Gene", all.x=FALSE, sort = FALSE)
head(cellmarkers)
dim(cellmarkers) #print all 41 genes assigned

bhybcountdatauqre_gene = merge(bhybcountdatauqre, ens_gene_names_chrpos_dedup_M20, by="id", all.x=FALSE, sort = FALSE)
head(bhybcountdatauqre_gene)
dim(bhybcountdatauqre_gene) 
bhybcountdatauqre_gene <- bhybcountdatauqre_gene[,c(16:19,1:15)]
write.table(bhybcountdatauqre_gene, "bhybcountdatauqre_gene.txt", sep = "\t", quote = F, append = F)

bhybcountdatauqre.norm_gene = merge(bhybcountdatauqre.norm, ens_gene_names_chrpos_dedup_M20, by="id", all.x=FALSE, sort = FALSE)
head(bhybcountdatauqre.norm_gene)
dim(bhybcountdatauqre.norm_gene) 
bhybcountdatauqre.norm_gene <- bhybcountdatauqre.norm_gene[,c(16:19,1:15)]
write.table(bhybcountdatauqre.norm_gene, "bhybcountdatauqre.norm_gene.txt", sep = "\t", quote = F, append = F)


bhybcountdatauqre.norm.lin = merge(cellmarkers,bhybcountdatauqre.norm, by="id", all.x=FALSE, sort = FALSE)
head(bhybcountdatauqre.norm.lin,2)
dim(bhybcountdatauqre.norm.lin) #Print all 41 genes overlapped and assigned

bhybcountdatauqre.norm.linre <- bhybcountdatauqre.norm.lin[,c(2, 7:20)]
head(bhybcountdatauqre.norm.linre,2)
rownames(bhybcountdatauqre.norm.linre) <- bhybcountdatauqre.norm.linre[,1]
head(bhybcountdatauqre.norm.linre)
bhybcountdatauqre.norm.linre1 <- bhybcountdatauqre.norm.linre[,c(2:15)]
head(bhybcountdatauqre.norm.linre1,2)
dim(bhybcountdatauqre.norm.linre1)
bhybcountdatauqre.norm.lineag <- as.matrix(bhybcountdatauqre.norm.linre1)



boxplot(log10(bhybcountdatauqre.norm.lineag), ylim=c(0,10), col = c("blue","blue","green","green","green","green","green","green","green","green","orange","orange","orange","orange"))
#save as bhybcountdatauqre.norm.lineag.sep.svg
library(EDASeq)
plotPCA(as.matrix(bhybcountdatauqre.norm.lineag))
plotPCA(as.matrix(bhybcountdatauqre.norm.lineag), labels=F, col =  c("cyan","red","blue","blue","darkgreen","darkgreen","blue","blue","darkgreen","darkgreen","navy","navy","red","red"))
ggsave("PCA_tbhybcountdatauqre.norm.lineag_scaleT.sep.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

write.table(bhybcountdatauqre.norm.lineag, "bhybcountdatauqre.norm.lineag.sep.txt", row.names = T, quote = F, append=F)
write.csv(bhybcountdatauqre.norm.lineag, "bhybcountdatauqre.norm.lineag.sep.csv")
#Dendrogram cluster
dim(bhybcountdatauqre.norm.lineag)  #41 14
# Compute distances and hierarchical clustering
bhybcountdatauqre.norm.lineag1 <- t(bhybcountdatauqre.norm.lineag[c(1:17,19:41),])
bhybcountdatauqre.norm.lineagscale <- scale(bhybcountdatauqre.norm.lineag1, center = TRUE, scale = TRUE)
ddbhybcountdatauqre.norm.lineagscale <- dist(as.matrix(bhybcountdatauqre.norm.lineagscale), method = "euclidean")
ddbhybcountdatauqre.norm.lineagscalehc <- hclust(ddbhybcountdatauqre.norm.lineagscale, method = "ward.D")
plot(ddbhybcountdatauqre.norm.lineagscalehc)
#save as ddbhybcountdatauqre.norm.lineagscalehc.sep.svg
head(bhybcountdatauqre.norm.lineag,1)                                                             

colfunc <- colorRampPalette(c("navy","white","firebrick3"))
head(bhybcountdatauqre.norm.lineag)
dim(bhybcountdatauqre.norm.lineag)  #41 14
summary(bhybcountdatauqre.norm.lineag)
library(pheatmap)
library(RColorBrewer)
breaksList2 = seq(0, 0.05, by = 0.001)
my_sample_row1 <- data.frame(Markers= c("ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Knockout"))
rownames(my_sample_row1) <- rownames(bhybcountdatauqre.norm.lineag)
my_colour1 = list(Markers = c(ESC = "#00ff00", Neural.cortex ="#008000", Endo.Meso = "#ff80b3", Knockout = "#B3A580"))
#pheatmap(selected_markers_T0_T12_1.ratio,treeheight_row = 0, cluster_cols=F, cluster_rows=T, treeheight_col = 0, gaps_col =NULL, gaps_row = NULL, border_color = "black", breaks = breaksList,color= colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)))
#Remove Otx2
pheatmap(bhybcountdatauqre.norm.lineag[c(1:17,19:41),],
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(length(breaksList2)),
         breaks = breaksList2,
         treeheight_row = 30,
         treeheight_col = 30,
         annotation_colors = my_colour1,
         fontsize = 8,
         annotation_row = my_sample_row1,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=T, 
         cluster_rows=F,
         cutree_cols = 3,
         cellwidth = 10, 
         cellheight = 10)
#save as filename = "pheatmap_bhybcountdatauqre.norm.lineag.sep.png")

dev.off()

#Z-score
z_Tbhybcountdatauqre.norm.lineag= scale(t(bhybcountdatauqre.norm.lineag), center = TRUE, scale = TRUE)
z_bhybcountdatauqre.norm.lineag <- t(z_Tbhybcountdatauqre.norm.lineag)
colfunc <- colorRampPalette(c("navy","white","firebrick3"))
head(z_bhybcountdatauqre.norm.lineag)
dim(z_bhybcountdatauqre.norm.lineag)  #41 14
summary(z_bhybcountdatauqre.norm.lineag)
library(pheatmap)
library(RColorBrewer)
breaksList1 = seq(-1, 1, by = 0.01)
my_sample_row1 <- data.frame(Markers= c("ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Knockout"))
rownames(my_sample_row1) <- rownames(z_bhybcountdatauqre.norm.lineag)
my_colour1 = list(Markers = c(ESC = "#ebeb14", Neural.cortex ="#8ca50e", Endo.Meso = "#ff80b3", Knockout = "#B3A580"))
#pheatmap(selected_markers_T0_T12_1.ratio,treeheight_row = 0, cluster_cols=F, cluster_rows=T, treeheight_col = 0, gaps_col =NULL, gaps_row = NULL, border_color = "black", breaks = breaksList,color= colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)))
#Remove Otx2
pheatmap(z_bhybcountdatauqre.norm.lineag[c(1:17,19:41),],
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(length(breaksList1)),
         breaks = breaksList1,
         treeheight_row = 30,
         treeheight_col = 30,
         annotation_colors = my_colour1,
         fontsize = 8,
         annotation_row = my_sample_row1,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=T, 
         cluster_rows=F,
         cutree_cols = 2,
         cellwidth = 10,
         cellheight = 10)

#save as filename = "pheatmap_z_bhybcountdatauqre.norm.lineag.sep.svg")

dev.off()
write.table(z_bhybcountdatauqre.norm.lineag, "z_bhybcountdatauqre.norm.lineag.txt", sep = "\t", append = F, quote = F)
#Select only specific genes (remove Otx2, Foxg1, Sox1, Mkx, Neurod1, Gfap and Pecam1.) as told by Basilia befroe z-score calculation
#Z-score
bhybcountdatauqre.norm.lineag_sel <- bhybcountdatauqre.norm.lineag[c(1:17,19,21:23,25,26,29:32,34:36,38:41),]
z_Tbhybcountdatauqre.norm.lineag_sel= scale(t(bhybcountdatauqre.norm.lineag_sel), center = TRUE, scale = TRUE)
z_bhybcountdatauqre.norm.lineag_sel <- t(z_Tbhybcountdatauqre.norm.lineag_sel)
colfunc <- colorRampPalette(c("navy","white","firebrick3"))
head(z_bhybcountdatauqre.norm.lineag_sel)
dim(z_bhybcountdatauqre.norm.lineag_sel)  #34 14
summary(z_bhybcountdatauqre.norm.lineag_sel)
library(pheatmap)
library(RColorBrewer)
breaksList1 = seq(-1, 1, by = 0.01)
my_sample_row1 <- data.frame(Markers= c("ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","ESC","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Neural.cortex","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Endo.Meso","Knockout"))
rownames(my_sample_row1) <- rownames(z_bhybcountdatauqre.norm.lineag_sel)
my_colour1 = list(Markers = c(ESC = "#ebeb14", Neural.cortex ="#8ca50e", Endo.Meso = "#ff80b3", Knockout = "#B3A580"))
#pheatmap(selected_markers_T0_T12_1.ratio,treeheight_row = 0, cluster_cols=F, cluster_rows=T, treeheight_col = 0, gaps_col =NULL, gaps_row = NULL, border_color = "black", breaks = breaksList,color= colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)))
#Remove Otx2
pheatmap(z_bhybcountdatauqre.norm.lineag_sel,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(length(breaksList1)),
         breaks = breaksList1,
         treeheight_row = 30,
         treeheight_col = 30,
         annotation_colors = my_colour1,
         fontsize = 8,
         annotation_row = my_sample_row1,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=T, 
         cluster_rows=F,
         cutree_cols = 2,
         cellwidth = 10,
         cellheight = 10)

#save as filename = "pheatmap_z_bhybcountdatauqre.norm.lineag_sel.sep.svg")

dev.off()

#No cluster
z_bhybcountdatauqre.norm.lineag_sel_re <- z_bhybcountdatauqre.norm.lineag_sel[,c(10,9,6,5,14:11,8,7,4:1)]
pheatmap(z_bhybcountdatauqre.norm.lineag_sel_re,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(length(breaksList1)),
         breaks = breaksList1,
         treeheight_row = 30,
         treeheight_col = 30,
         annotation_colors = my_colour1,
         fontsize = 8,
         annotation_row = my_sample_row1,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=F, 
         cluster_rows=F,
         cutree_cols = 2,
         cellwidth = 10,
         cellheight = 10)

#save as filename = "pheatmap_z_bhybcountdatauqre.norm.lineag_sel_re.sep.svg")

#Export files
dim(bhybcountdatauqre.norm.lineag_sel)
bhybcountdatauqre.norm.lineag_sel_gene <- data.frame(bhybcountdatauqre.norm.lineag_sel)
bhybcountdatauqre.norm.lineag_sel_gene$Gene <- rownames(bhybcountdatauqre.norm.lineag_sel_gene)
bhybcountdatauqre.norm.lineag_sel_gene <- bhybcountdatauqre.norm.lineag_sel_gene[,c(15,1:14)]
head(bhybcountdatauqre.norm.lineag_sel_gene)

dim(z_bhybcountdatauqre.norm.lineag_sel)
z_bhybcountdatauqre.norm.lineag_sel_gene <- data.frame(z_bhybcountdatauqre.norm.lineag_sel)
z_bhybcountdatauqre.norm.lineag_sel_gene$Gene <- rownames(z_bhybcountdatauqre.norm.lineag_sel_gene)
z_bhybcountdatauqre.norm.lineag_sel_gene <- z_bhybcountdatauqre.norm.lineag_sel_gene[,c(15,1:14)]
head(z_bhybcountdatauqre.norm.lineag_sel_gene)

combinedbhybcountdatauqre.norm.lineag_sel_gene <-  cbind.data.frame(bhybcountdatauqre.norm.lineag_sel_gene,
                                                                    z_bhybcountdatauqre.norm.lineag_sel_gene)


head(combinedbhybcountdatauqre.norm.lineag_sel_gene)
library(WriteXLS)
library("writexl")
write_xlsx(combinedbhybcountdatauqre.norm.lineag_sel_gene, "combinedbhybcountdatauqre.norm.lineag_sel_gene.xlsx")
###############################---------------------------------------------------
#Extract imprinted genes coordinates
#created after basila imprinted gene selection /media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/selected.markers
imprinted.gene.list <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb/imprinted_gene_name.txt", header = F)
head(imprinted.gene.list)
colnames(imprinted.gene.list) <- c("Gene")
dim(imprinted.gene.list) #print 188 genes

imprinted.genes = merge(imprinted.gene.list, ens_gene_names_chrpos_dedup_M20, by="Gene", all.x=FALSE, sort = FALSE)
head(imprinted.genes)
dim(imprinted.genes) #print all 158 genes assigned


bhybcountdatauqre.norm.imp = merge(imprinted.genes,bhybcountdatauqre.norm, by="id", all.x=FALSE, sort = FALSE)
head(bhybcountdatauqre.norm.imp,2)
dim(bhybcountdatauqre.norm.imp) #Print all 158 genes overlapped and assigned

#Select samples and gene symbol
bhybcountdatauqre.norm.impre <- bhybcountdatauqre.norm.imp[,c(2, 6:7,16:19)]
head(bhybcountdatauqre.norm.impre,2)
rownames(bhybcountdatauqre.norm.impre) <- bhybcountdatauqre.norm.impre[,1]
head(bhybcountdatauqre.norm.impre,2)
bhybcountdatauqre.norm.impre1 <- bhybcountdatauqre.norm.impre[,c(2:7)]
head(bhybcountdatauqre.norm.impre1,2)
dim(bhybcountdatauqre.norm.impre1)
bhybcountdatauqre.norm.imprint <- as.matrix(bhybcountdatauqre.norm.impre1)



boxplot(log10(bhybcountdatauqre.norm.imprint), ylim=c(0,10), col = c("darkgreen","darkred","green","green","red","red"))
#save as bhybcountdatauqre.norm.imprint.sep.svg
library(EDASeq)
plotPCA(as.matrix(bhybcountdatauqre.norm.imprint))
plotPCA(as.matrix(bhybcountdatauqre.norm.imprint), labels=F, col =  c("darkgreen","darkred","green","green","red","red"))
ggsave("PCA_tbhybcountdatauqre.norm.imprint_scaleT.sep.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

write.table(bhybcountdatauqre.norm.imprint, "bhybcountdatauqre.norm.imprint.sep.txt", row.names = T, quote = F, append=F)
#Dendrogram cluster
dim(bhybcountdatauqre.norm.imprint)  #158 6

# Compute distances and hierarchical clustering
bhybcountdatauqre.norm.imprintT <- t(bhybcountdatauqre.norm.imprint)
bhybcountdatauqre.norm.imprintscale <- scale(bhybcountdatauqre.norm.imprintT, center = TRUE, scale = TRUE)
ddbhybcountdatauqre.norm.imprintscale <- dist(as.matrix(bhybcountdatauqre.norm.imprintscale), method = "euclidean")
ddbhybcountdatauqre.norm.imprintscalehc <- hclust(ddbhybcountdatauqre.norm.imprintscale, method = "ward.D")
plot(ddbhybcountdatauqre.norm.imprintscalehc)
#save as ddbhybcountdatauqre.norm.imprintscalehc.sep.svg
head(bhybcountdatauqre.norm.imprint,1)                                                             

#save as ddbhybcountdatauqre.norm.imprintscalehc.sep.svg
head(bhybcountdatauqre.norm.imprint,1)                                                             

colfunc <- colorRampPalette(c("navy","white","firebrick3"))
head(bhybcountdatauqre.norm.imprint)
dim(bhybcountdatauqre.norm.imprint)  #158 6
summary(bhybcountdatauqre.norm.imprint)
library(pheatmap)
library(RColorBrewer)
breaksList2 = seq(0, 0.05, by = 0.001)
my_imp_colour1 = list(Markers = c(Imp = "#00ff00"))

pheatmap(bhybcountdatauqre.norm.imprint,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(length(breaksList2)),
         breaks = breaksList2,
         treeheight_row = 30,
         treeheight_col = 30,
         annotation_colors = my_imp_colour1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=T, 
         cluster_rows=F,
         cutree_cols = 3,
         cellwidth = 10, 
         cellheight = 10)
#save as filename = "pheatmap_bhybcountdatauqre.norm.imprint.sep.svg")

dev.off()

#Z-score
z_Tbhybcountdatauqre.norm.imprint= scale(t(bhybcountdatauqre.norm.imprint), center = TRUE, scale = TRUE)
z_bhybcountdatauqre.norm.imprint <- t(z_Tbhybcountdatauqre.norm.imprint)
colfunc <- colorRampPalette(c("navy","white","firebrick3"))
head(z_bhybcountdatauqre.norm.imprint)
dim(z_bhybcountdatauqre.norm.imprint)  #158 6
summary(z_bhybcountdatauqre.norm.imprint)
library(pheatmap)
library(RColorBrewer)
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_bhybcountdatauqre.norm.imprint,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(length(breaksList1)),
         breaks = breaksList1,
         treeheight_row = 30,
         treeheight_col = 30,
         annotation_colors = my_imp_colour1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=T, 
         cluster_rows=F,
         cutree_cols = 2,
         cellwidth = 10,
         cellheight = 10)

#save as filename = "pheatmap_z_bhybcountdatauqre.norm.imprint.sep.svg")


#or
pheatmap(z_bhybcountdatauqre.norm.imprint,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(length(breaksList1)),
         breaks = breaksList1,
         treeheight_row = 30,
         treeheight_col = 30,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=T, 
         cluster_rows=F,
         cutree_cols = 2,
         cellwidth = 10,
         cellheight = 10)
#save as filename = "pheatmap_z_bhybcountdatauqre.norm.imprintall.sep.svg")

#Dendrogram of all genes

#Dendrogram cluster
dim(bhybcountdatauqre.norm)
head(bhybcountdatauqre.norm)
# Compute distances and hierarchical clustering for only our sample
bhybcountdatauqre.normT <- t(bhybcountdatauqre.norm[,c(2,3,12:15)])
bhybcountdatauqre.normscale <- scale(bhybcountdatauqre.normT, center = TRUE, scale = TRUE)
ddbhybcountdatauqre.normscale <- dist(as.matrix(bhybcountdatauqre.normscale), method = "euclidean")
ddbhybcountdatauqre.normscalehc <- hclust(ddbhybcountdatauqre.normscale, method = "ward.D")
plot(ddbhybcountdatauqre.normscalehc)
#save as ddbhybcountdatauqre.normscalehc.sep.svg

dev.off()

#------------------
#Extract limited imprinted genes coordinates
#created after basila imprinted gene selection /media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/selected.markers
#imprinted.gene.limited <- read.table("imprinted_gene_name_limited.txt", header = F)
imprinted.gene.limited <- read.table("imprinted_gene_name_limited.txt", header = F)

head(imprinted.gene.limited)
colnames(imprinted.gene.limited) <- c("Gene")
dim(imprinted.gene.limited) #print 39 genes

imprinted.genes_lim_chr = merge(imprinted.gene.limited, ens_gene_names_chrpos_dedup_M20, by="Gene", all.x=FALSE, sort = FALSE)
head(imprinted.genes_lim_chr)
dim(imprinted.genes_lim_chr) #print all 39 genes assigned


bhybcountdatauqre.norm_limited.imp = merge(imprinted.genes_lim_chr,bhybcountdatauqre.norm, by="id", all.x=FALSE, sort = FALSE)
head(bhybcountdatauqre.norm_limited.imp,2)
dim(bhybcountdatauqre.norm_limited.imp) #Print all 39 genes overlapped and assigned

bhybcountdatauqre.norm_limited.impre <- bhybcountdatauqre.norm_limited.imp[,c(2, 6:7,16:19)]
head(bhybcountdatauqre.norm_limited.impre,2)
rownames(bhybcountdatauqre.norm_limited.impre) <- bhybcountdatauqre.norm_limited.impre[,1]
head(bhybcountdatauqre.norm_limited.impre,2)
bhybcountdatauqre.norm_limited.impre1 <- bhybcountdatauqre.norm_limited.impre[,c(2:7)]
head(bhybcountdatauqre.norm_limited.impre1,2)
dim(bhybcountdatauqre.norm_limited.impre1)
bhybcountdatauqre.norm_limited.imprint <- as.matrix(bhybcountdatauqre.norm_limited.impre1)

boxplot(log10(bhybcountdatauqre.norm_limited.imprint), ylim=c(0,10), col = c("darkgreen","darkred","green","green","red","red"))
#save as bhybcountdatauqre.norm_limited.imprint.sep.svg
library(EDASeq)
plotPCA(as.matrix(bhybcountdatauqre.norm_limited.imprint))
plotPCA(as.matrix(bhybcountdatauqre.norm_limited.imprint), labels=F, col =  c("darkgreen","darkred","green","green","red","red"))
ggsave("PCA_tbhybcountdatauqre.norm_limited.imprint_scaleT.sep.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

write.table(bhybcountdatauqre.norm_limited.imprint, "bhybcountdatauqre.norm_limited.imprint.sep.txt", row.names = T, quote = F, append=F)
#Dendrogram cluster
dim(bhybcountdatauqre.norm_limited.imprint)  #39 6

# Compute distances and hierarchical clustering
bhybcountdatauqre.norm_limited.imprintT <- t(bhybcountdatauqre.norm_limited.imprint)
bhybcountdatauqre.norm_limited.imprintscale <- scale(bhybcountdatauqre.norm_limited.imprintT, center = TRUE, scale = TRUE)
ddbhybcountdatauqre.norm_limited.imprintscale <- dist(as.matrix(bhybcountdatauqre.norm_limited.imprintscale), method = "euclidean")
ddbhybcountdatauqre.norm_limited.imprintscalehc <- hclust(ddbhybcountdatauqre.norm_limited.imprintscale, method = "ward.D")
plot(ddbhybcountdatauqre.norm_limited.imprintscalehc)
#save as ddbhybcountdatauqre.norm_limited.imprintscalehc.sep.svg
head(bhybcountdatauqre.norm_limited.imprint,1)                                                             

#save as ddbhybcountdatauqre.norm_limited.imprintscalehc.sep.svg
head(bhybcountdatauqre.norm_limited.imprint,1)                                                             

colfunc <- colorRampPalette(c("navy","white","firebrick3"))
head(bhybcountdatauqre.norm_limited.imprint)
dim(bhybcountdatauqre.norm_limited.imprint)  #158 6
summary(bhybcountdatauqre.norm_limited.imprint)
library(pheatmap)
library(RColorBrewer)
breaksList2 = seq(0, 0.05, by = 0.001)
my_imp_colour1 = list(Markers = c(Imp = "#00ff00"))

pheatmap(bhybcountdatauqre.norm_limited.imprint,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(length(breaksList2)),
         breaks = breaksList2,
         treeheight_row = 30,
         treeheight_col = 30,
         annotation_colors = my_imp_colour1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=T, 
         cluster_rows=F,
         cutree_cols = 3,
         cellwidth = 10, 
         cellheight = 10)
#save as filename = "pheatmap_bhybcountdatauqre.norm_limited.imprint.sep.svg")

dev.off()

#Z-score
z_Tbhybcountdatauqre.norm_limited.imprint= scale(t(bhybcountdatauqre.norm_limited.imprint), center = TRUE, scale = TRUE)
z_bhybcountdatauqre.norm_limited.imprint <- t(z_Tbhybcountdatauqre.norm_limited.imprint)
colfunc <- colorRampPalette(c("navy","white","firebrick3"))
head(z_bhybcountdatauqre.norm_limited.imprint)
dim(z_bhybcountdatauqre.norm_limited.imprint)  #158 6
summary(z_bhybcountdatauqre.norm_limited.imprint)
library(pheatmap)
library(RColorBrewer)
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_bhybcountdatauqre.norm_limited.imprint,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(length(breaksList1)),
         breaks = breaksList1,
         treeheight_row = 30,
         treeheight_col = 30,
         annotation_colors = my_imp_colour1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=T, 
         cluster_rows=F,
         cutree_cols = 2,
         cellwidth = 10,
         cellheight = 10)

#save as filename = "pheatmap_z_bhybcountdatauqre.norm_limited.imprint.sep.svg")



#Z-score
dim(bhybcountdatauqre.norm)
head(bhybcountdatauqre.norm)
z_Tbhybcountdatauqre.norm= scale(t(bhybcountdatauqre.norm[,2:15]), center = TRUE, scale = TRUE)
z_bhybcountdatauqre.norm <- t(z_Tbhybcountdatauqre.norm)
head(z_bhybcountdatauqre.norm)
dim(z_bhybcountdatauqre.norm)  #54750 14
summary(z_bhybcountdatauqre.norm)
z_bhybcountdatauqre.norm <- data.frame(z_bhybcountdatauqre.norm)
z_bhybcountdatauqre.norm["id"] <- rownames(z_bhybcountdatauqre.norm)
z_bhybcountdatauqre.norm_gene = merge(z_bhybcountdatauqre.norm, ens_gene_names_chrpos_dedup_M20, by="id", all.x=FALSE, sort = FALSE)
head(z_bhybcountdatauqre.norm_gene)
dim(z_bhybcountdatauqre.norm_gene) 
z_bhybcountdatauqre.norm_gene <- z_bhybcountdatauqre.norm_gene[,c(16:19,1:15)]
write.table(z_bhybcountdatauqre.norm_gene, "z_bhybcountdatauqre.norm_gene.txt", sep = "\t", quote = F, append = F)



head(bhybcountdatauqre.norm[,c(2:3,12:15)],2)
dim(bhybcountdatauqre.norm[,c(2:3,12:15)])

#Z-score
z_Tbhybcountdatauqre.int_norm = scale(t(bhybcountdatauqre.norm[,c(2:3,12:15)]), center = TRUE, scale = TRUE)
z_bhybcountdatauqre.int_norm <- t(z_Tbhybcountdatauqre.int_norm)
head(z_bhybcountdatauqre.int_norm)
dim(z_bhybcountdatauqre.int_norm)  #54750  6
summary(z_bhybcountdatauqre.int_norm)
z_bhybcountdatauqre.int_norm <- data.frame(z_bhybcountdatauqre.int_norm)
z_bhybcountdatauqre.int_norm["id"] <- rownames(z_bhybcountdatauqre.int_norm)
z_bhybcountdatauqre.int_norm_gene = merge(z_bhybcountdatauqre.int_norm, ens_gene_names_chrpos_dedup_M20, by="id", all.x=FALSE, sort = FALSE)
head(z_bhybcountdatauqre.int_norm_gene)
dim(z_bhybcountdatauqre.int_norm_gene) 
z_bhybcountdatauqre.int_norm_gene <- z_bhybcountdatauqre.int_norm_gene[,c(8:11,1:7)]
write.table(z_bhybcountdatauqre.int_norm_gene, "z_bhybcountdatauqre.int_norm_gene.txt", sep = "\t", quote = F, append = F, row.names = F)


#------------------
#Extract limited imprrinted genes coordinates
#created after basila imprrinted gene selection /media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/selected.markers
#imprrinted.gene.limited <- read.table("imprrinted_gene_name_limited.txt", header = F)
imprrinted.gene.limited <- read.table("imprinted_gene_name_limited2.txt", header = F)

head(imprrinted.gene.limited)
colnames(imprrinted.gene.limited) <- c("Gene")
dim(imprrinted.gene.limited) #print 40 genes

imprrinted.genes_lim_chr = merge(imprrinted.gene.limited, ens_gene_names_chrpos_dedup_M20, by="Gene", all.x=FALSE, sort = FALSE)
head(imprrinted.genes_lim_chr)
dim(imprrinted.genes_lim_chr) #print all 40 genes assigned


bhybcountdatauqre.norm_limited.impr = merge(imprrinted.genes_lim_chr,bhybcountdatauqre.norm, by="id", all.x=FALSE, sort = FALSE)
head(bhybcountdatauqre.norm_limited.impr,2)
dim(bhybcountdatauqre.norm_limited.impr) #Print all 40 genes overlapped and assigned

bhybcountdatauqre.norm_limited.imprre <- bhybcountdatauqre.norm_limited.impr[,c(2, 6:19)]
head(bhybcountdatauqre.norm_limited.imprre,2)
rownames(bhybcountdatauqre.norm_limited.imprre) <- bhybcountdatauqre.norm_limited.imprre[,1]
head(bhybcountdatauqre.norm_limited.imprre,2)
bhybcountdatauqre.norm_limited.imprre1 <- bhybcountdatauqre.norm_limited.imprre[,2:15]
head(bhybcountdatauqre.norm_limited.imprre1,2)
dim(bhybcountdatauqre.norm_limited.imprre1)
bhybcountdatauqre.norm_limited.imprrint <- as.matrix(bhybcountdatauqre.norm_limited.imprre1)

boxplot(log10(bhybcountdatauqre.norm_limited.imprrint), ylim=c(0,10))
#save as bhybcountdatauqre.norm_limited.imprrint.sep.svg
library(EDASeq)
plotPCA(as.matrix(bhybcountdatauqre.norm_limited.imprrint))
ggsave("PCA_tbhybcountdatauqre.norm_limited.imprrint_scaleT.sep.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

write.table(bhybcountdatauqre.norm_limited.imprrint, "bhybcountdatauqre.norm_limited.imprrint.sep.txt", row.names = T, quote = F, append=F)
#Dendrogram cluster
dim(bhybcountdatauqre.norm_limited.imprrint)  #40 14

# Compute distances and hierarchical clustering
bhybcountdatauqre.norm_limited.imprrintT <- t(bhybcountdatauqre.norm_limited.imprrint)
bhybcountdatauqre.norm_limited.imprrintscale <- scale(bhybcountdatauqre.norm_limited.imprrintT, center = TRUE, scale = TRUE)
ddbhybcountdatauqre.norm_limited.imprrintscale <- dist(as.matrix(bhybcountdatauqre.norm_limited.imprrintscale), method = "euclidean")
ddbhybcountdatauqre.norm_limited.imprrintscalehc <- hclust(ddbhybcountdatauqre.norm_limited.imprrintscale, method = "ward.D")
plot(ddbhybcountdatauqre.norm_limited.imprrintscalehc)
#save as ddbhybcountdatauqre.norm_limited.imprrintscalehc.sep.svg
head(bhybcountdatauqre.norm_limited.imprrint,1)                                                             

#save as ddbhybcountdatauqre.norm_limited.imprrintscalehc.sep.svg
head(bhybcountdatauqre.norm_limited.imprrint,1)                                                             

colfunc <- colorRampPalette(c("navy","white","firebrick3"))
head(bhybcountdatauqre.norm_limited.imprrint)
dim(bhybcountdatauqre.norm_limited.imprrint)  #40  14
summary(bhybcountdatauqre.norm_limited.imprrint)
library(pheatmap)
library(RColorBrewer)
breaksList2 = seq(0, 0.05, by = 0.001)
my_impr_colour1 = list(Markers = c(impr = "#00ff00"))

pheatmap(bhybcountdatauqre.norm_limited.imprrint,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(length(breaksList2)),
         breaks = breaksList2,
         treeheight_row = 30,
         treeheight_col = 30,
         annotation_colors = my_impr_colour1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=T, 
         cluster_rows=F,
         cutree_cols = 3,
         cellwidth = 10, 
         cellheight = 10)
#save as filename = "pheatmap_bhybcountdatauqre.norm_limited.imprrint.sep.svg")

dev.off()

#Z-score
z_Tbhybcountdatauqre.norm_limited.imprrint= scale(t(bhybcountdatauqre.norm_limited.imprrint), center = TRUE, scale = TRUE)
z_bhybcountdatauqre.norm_limited.imprrint <- t(z_Tbhybcountdatauqre.norm_limited.imprrint)
colfunc <- colorRampPalette(c("navy","white","firebrick3"))
head(z_bhybcountdatauqre.norm_limited.imprrint)
dim(z_bhybcountdatauqre.norm_limited.imprrint)  #40  14
summary(z_bhybcountdatauqre.norm_limited.imprrint)
library(pheatmap)
library(RColorBrewer)
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_bhybcountdatauqre.norm_limited.imprrint,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(length(breaksList1)),
         breaks = breaksList1,
         treeheight_row = 30,
         treeheight_col = 30,
         annotation_colors = my_impr_colour1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=T, 
         cluster_rows=F,
         cutree_cols = 2,
         cellwidth = 10,
         cellheight = 10)

#save as filename = "pheatmap_z_bhybcountdatauqre.norm_limited.imprrint.sep.svg")


head(bhybcountdatauqre.norm_limited.imprrint[,c(1,2,11:14)])
#Z-score
z_Tbhybcountdatauqre.int_norm_limited.imprrint= scale(t(bhybcountdatauqre.norm_limited.imprrint[,c(1,2,11:14)]), center = TRUE, scale = TRUE)
z_bhybcountdatauqre.int_norm_limited.imprrint <- t(z_Tbhybcountdatauqre.int_norm_limited.imprrint)
colfunc <- colorRampPalette(c("navy","white","firebrick3"))
head(z_bhybcountdatauqre.int_norm_limited.imprrint)
dim(z_bhybcountdatauqre.int_norm_limited.imprrint)  #40  6
summary(z_bhybcountdatauqre.int_norm_limited.imprrint)
library(pheatmap)
library(RColorBrewer)
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_bhybcountdatauqre.int_norm_limited.imprrint,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(length(breaksList1)),
         breaks = breaksList1,
         treeheight_row = 30,
         treeheight_col = 30,
         annotation_colors = my_impr_colour1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "white",
         cluster_cols=T, 
         cluster_rows=F,
         cutree_cols = 2,
         cellwidth = 10,
         cellheight = 10)

#save as filename = "pheatmap_z_bhybcountdatauqre.int_norm_limited.imprrint.sep.svg")

write.table(bhybcountdatauqre.norm_limited.imprrint, "bhybcountdatauqre.norm_limited.imprrint.txt", quote=F, append = F, row.names = T, sep="\t")

write.table(bhybcountdata, "bhybcountdata.txt", sep="\t", quote = FALSE, append = FALSE)

forgenelength <- read.table("Bulk_hyb_ourT0j_star-featureCounts_GRCm38.mm10.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
forgenelength <- forgenelength[ ,5:ncol(forgenelength)]
head(forgenelength)
gene.length <- forgenelength$Length

head(bhybcountdata)
countforfpkm  <-  bhybcountdata
head(countforfpkm)
countforfpkm <- countforfpkm[order(rownames(countforfpkm)),]



countforfpkm_filt1 <- data.frame(countforfpkm)
head(countforfpkm_filt1)

#NOISeq
myfactors = read.table("myfactors.txt", header=TRUE)
head(myfactors)
rownames(myfactors)
colnames(myfactors)

#Plot first PCA without filtering or normalisation of data: see ggplot2 
mydata1 <- NOISeq::readData(data=countforfpkm_filt1, length = gene.length, factors=myfactors)
myPCA1 = dat(mydata1, type = "PCA")
par(mfrow = c(1, 2))
explo.plot(myPCA1, factor = "btypes")

#Normalisation
myRPKM = NOISeq::rpkm(assayData(mydata1)$exprs, long = gene.length, k = 0, lc = 1)
head(myRPKM)
mydata3 <- NOISeq::readData(data=myRPKM, factors=myfactors)
myPCA3 = dat(mydata3, type = "PCA")
par(mfrow = c(1, 2))
explo.plot(myPCA3, factor = "btypes")

#Low-count filtering
dim(myRPKM)
myRPKMFilt <- myRPKM[rowSums(myRPKM) > 0, ]
dim(myRPKMFilt)
myRPKMFilt_sorted = myRPKMFilt[order(rownames(myRPKMFilt)),]
head(myRPKMFilt_sorted)
dim(myRPKMFilt_sorted)
write.table(myRPKMFilt_sorted, "myRPKMFilt_sorted.txt", sep="\t", quote = FALSE, append = FALSE)



#edgeR rpkm()
library(edgeR)
dgeforfpkm <- DGEList(counts=countforfpkm_filt1,genes=data.frame(Length=gene.length), group = c("WT","ZFP57KO","WT","WT","WT","WT","WT","WT","WT","WT","WT","WT","ZFP57KO","ZFP57KO"))
dgeforfpkm <- calcNormFactors(dgeforfpkm)
dgeFPKM <- rpkm(dgeforfpkm, dgeforfpkm$genes$Length)
head(dgeFPKM)
write.table(dgeFPKM, "dgeFPKM.txt", sep="\t", quote = FALSE, append = FALSE)
dgeFPKM <- data.frame(dgeFPKM)
dgeFPKM["id"] <- rownames(dgeFPKM)
dgeFPKM_gene = merge(dgeFPKM, ens_gene_names_chrpos_dedup_M20, by="id", all.x=FALSE, sort = FALSE)
head(dgeFPKM_gene)
dim(dgeFPKM_gene)

#Check rpkm conversion
dim(countforfpkm_filt1)
sum(countforfpkm_filt1$ourT0WTJB1) #41435450
head(countforfpkm_filt1,1) #8282
head(dgeFPKM_gene,1) #191.5156
head(dgeforfpkm$genes$Length) #1070
#Checked: ((readcounts * 10^6 * 10^3)/(lib.size * genelength))/normfactor = ((8282 * 1000000000)/(41435450 * 1070))/0.9753831 = 191.5256

#Take average of T12 wt n Ko
head(dgeFPKM_gene)
dgeFPKM_gene["ourT12WTavg"] <- (dgeFPKM_gene$ourT12WT1 + dgeFPKM_gene$ourT12WT2)/2
dgeFPKM_gene["ourT12ZFP57KOavg"] <- (dgeFPKM_gene$ourT12ZFP57KO1 + dgeFPKM_gene$ourT12ZFP57KO2)/2

ourT0WTJB1_FPKM <- dgeFPKM_gene[,c(2,19)]
ourT0WTJB1_FPKM <- ourT0WTJB1_FPKM[order(-ourT0WTJB1_FPKM$ourT0WTJB1),]
ourT0WTJB1_FPKM["Rank"] <- paste0(seq(1:length(ourT0WTJB1_FPKM$Gene)))
ourT0WTJB1_FPKM["SampleType"] <- "ourT0WTJB1"
head(ourT0WTJB1_FPKM)
write.table(ourT0WTJB1_FPKM, "ourT0WTJB1_FPKM.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)


ourT0ZFP57KOJB1_FPKM <- dgeFPKM_gene[,c(3,19)]
ourT0ZFP57KOJB1_FPKM <- ourT0ZFP57KOJB1_FPKM[order(-ourT0ZFP57KOJB1_FPKM$ourT0ZFP57KOJB1),]
ourT0ZFP57KOJB1_FPKM["Rank"] <- paste0(seq(1:length(ourT0ZFP57KOJB1_FPKM$Gene)))
ourT0ZFP57KOJB1_FPKM["SampleType"] <- "ourT0ZFP57KOJB1"
head(ourT0ZFP57KOJB1_FPKM)
write.table(ourT0ZFP57KOJB1_FPKM, "ourT0ZFP57KOJB1_FPKM.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)

ourT12WT1_FPKM <- dgeFPKM_gene[,c(12,19)]
ourT12WT1_FPKM <- ourT12WT1_FPKM[order(-ourT12WT1_FPKM$ourT12WT1),]
ourT12WT1_FPKM["Rank"] <- paste0(seq(1:length(ourT12WT1_FPKM$Gene)))
ourT12WT1_FPKM["SampleType"] <- "ourT12WT1"
head(ourT12WT1_FPKM)
write.table(ourT12WT1_FPKM, "ourT12WT1_FPKM.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)


ourT12WT2_FPKM <- dgeFPKM_gene[,c(13,19)]
ourT12WT2_FPKM <- ourT12WT2_FPKM[order(-ourT12WT2_FPKM$ourT12WT2),]
ourT12WT2_FPKM["Rank"] <- paste0(seq(1:length(ourT12WT2_FPKM$Gene)))
ourT12WT2_FPKM["SampleType"] <- "ourT12WT2"
head(ourT12WT2_FPKM)
write.table(ourT12WT2_FPKM, "ourT12WT2_FPKM.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)

ourT12ZFP57KO1_FPKM <- dgeFPKM_gene[,c(14,19)]
ourT12ZFP57KO1_FPKM <- ourT12ZFP57KO1_FPKM[order(-ourT12ZFP57KO1_FPKM$ourT12ZFP57KO1),]
ourT12ZFP57KO1_FPKM["Rank"] <- paste0(seq(1:length(ourT12ZFP57KO1_FPKM$Gene)))
ourT12ZFP57KO1_FPKM["SampleType"] <- "ourT12ZFP57KO1"
head(ourT12ZFP57KO1_FPKM)
write.table(ourT12ZFP57KO1_FPKM, "ourT12ZFP57KO1_FPKM.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)


ourT12ZFP57KO2_FPKM <- dgeFPKM_gene[,c(15,19)]
ourT12ZFP57KO2_FPKM <- ourT12ZFP57KO2_FPKM[order(-ourT12ZFP57KO2_FPKM$ourT12ZFP57KO2),]
ourT12ZFP57KO2_FPKM["Rank"] <- paste0(seq(1:length(ourT12ZFP57KO2_FPKM$Gene)))
ourT12ZFP57KO2_FPKM["SampleType"] <- "ourT12ZFP57KO2"
head(ourT12ZFP57KO2_FPKM)
write.table(ourT12ZFP57KO2_FPKM, "ourT12ZFP57KO2_FPKM.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)
dim(dgeFPKM_gene)
ourT12WTavg_FPKM <- dgeFPKM_gene[,c(20,19)]
ourT12WTavg_FPKM <- ourT12WTavg_FPKM[order(-ourT12WTavg_FPKM$ourT12WTavg),]
ourT12WTavg_FPKM["Rank"] <- paste0(seq(1:length(ourT12WTavg_FPKM$Gene)))
ourT12WTavg_FPKM["SampleType"] <- "ourT12WTavg"
head(ourT12WTavg_FPKM)
write.table(ourT12WTavg_FPKM, "ourT12WTavg_FPKM.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)

ourT12ZFP57KOavg_FPKM <- dgeFPKM_gene[,c(21,19)]
ourT12ZFP57KOavg_FPKM <- ourT12ZFP57KOavg_FPKM[order(-ourT12ZFP57KOavg_FPKM$ourT12ZFP57KOavg),]
ourT12ZFP57KOavg_FPKM["Rank"] <- paste0(seq(1:length(ourT12ZFP57KOavg_FPKM$Gene)))
ourT12ZFP57KOavg_FPKM["SampleType"] <- "ourT12ZFP57KOavg"
head(ourT12ZFP57KOavg_FPKM)
write.table(ourT12ZFP57KOavg_FPKM, "ourT12ZFP57KOavg_FPKM.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F)

fgrep -f imprinted_gene_name_limited2.txt ourT0WTJB1_FPKM.txt -w | grep "mt-Th" -v | sort -k2,2 > ourT0WTJB1_FPKM_Imp_lim.txt
fgrep -f imprinted_gene_name_limited2.txt ourT0ZFP57KOJB1_FPKM.txt -w | grep "mt-Th" -v | sort -k2,2 > ourT0ZFP57KOJB1_FPKM_Imp_lim.txt
fgrep -f imprinted_gene_name_limited2.txt ourT12WT1_FPKM.txt -w | grep "mt-Th" -v | sort -k2,2 > ourT12WT1_FPKM_Imp_lim.txt
fgrep -f imprinted_gene_name_limited2.txt ourT12WT2_FPKM.txt -w | grep "mt-Th" -v | sort -k2,2 > ourT12WT2_FPKM_Imp_lim.txt
fgrep -f imprinted_gene_name_limited2.txt ourT12ZFP57KO1_FPKM.txt -w | grep "mt-Th" -v | sort -k2,2 > ourT12ZFP57KO1_FPKM_Imp_lim.txt
fgrep -f imprinted_gene_name_limited2.txt ourT12ZFP57KO2_FPKM.txt -w | grep "mt-Th" -v | sort -k2,2 > ourT12ZFP57KO2_FPKM_Imp_lim.txt

fgrep -f imprinted_gene_name_limited2.txt ourT12WTavg_FPKM.txt -w | grep "mt-Th" -v | sort -k2,2 > ourT12WTavg_FPKM_Imp_lim.txt
fgrep -f imprinted_gene_name_limited2.txt ourT12ZFP57KOavg_FPKM.txt -w | grep "mt-Th" -v | sort -k2,2 > ourT12ZFP57KOavg_FPKM_Imp_lim.txt

paste *Imp_lim.txt  > t0_t12_rank_FPKM_Imp_lim.txt

fgrep -f imprinted_gene_name_limited2.txt  o*FPKM.txt  -w   > imprinted_gene_name_o_FPKM.txt
t0_t12_rank_FPKM_Imp_lim <- read.table("t0_t12_rank_FPKM_Imp_lim.txt", header = F, stringsAsFactors = F)
rownames(t0_t12_rank_FPKM_Imp_lim) <- t0_t12_rank_FPKM_Imp_lim$V2
head(t0_t12_rank_FPKM_Imp_lim,2)
tail(t0_t12_rank_FPKM_Imp_lim,2)

t0_t12_rank_FPKM_Imp_lim_re <- t0_t12_rank_FPKM_Imp_lim[,c(3,7,11,15,19,23,27,31)]
head(t0_t12_rank_FPKM_Imp_lim_re)
colnames(t0_t12_rank_FPKM_Imp_lim_re) <- c("ourT0WTJB1","ourT0ZFP57KOJB1","ourT12WT1","ourT12WT2","ourT12WTavg","ourT12ZFP57KO1","ourT12ZFP57KO2","ourT12ZFP57KOavg")
head(t0_t12_rank_FPKM_Imp_lim_re)
dim(t0_t12_rank_FPKM_Imp_lim_re)
write.table(t0_t12_rank_FPKM_Imp_lim_re, "t0_t12_rank_FPKM_Imp_lim_re.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)
with(t0_t12_rank_FPKM_Imp_lim_re, plot(ourT0WTJB1, ourT0ZFP57KOJB1, pch=21, cex = 0.2, main="Scatter Plot", xlim=c(0,54750), ylim=c(0,54750),lwd = 0.3,col="#525254", bg="grey",ylab="JB1 Zfp57-/- Rank", xlab="JB1 WT Rank", bty = 'n'))
with(subset(t0_t12_rank_FPKM_Imp_lim_re), points(ourT0WTJB1, ourT0ZFP57KOJB1, pch=21, cex = 0.8, xlim=c(0,54750), ylim=c(0,54750), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') #reddish
with(subset(t0_t12_rank_FPKM_Imp_lim_re), points(ourT12WTavg, ourT12ZFP57KOavg, pch=21, cex = 0.8, xlim=c(0,54750), ylim=c(0,54750), lwd = 0.5,col="red", bg="red"),ylab="", xlab="", bty = 'n') #reddish

lines(x = c(0,54750), y = c(0,54750), col = "#464647", lty = 2)
#Save manually as t0_t12_rank_FPKM_Imp_lim_re.svg


plot(t0_t12_rank_FPKM_Imp_lim_re$ourT0WTJB1, col="darkgreen", pch =16, cex = 0.8, xlim=c(0,40), ylim=c(0,54750))
par(new=T)
plot(t0_t12_rank_FPKM_Imp_lim_re$ourT0ZFP57KOJB1, col="darkred", pch =16, cex = 0.8, xlim=c(0,40), ylim=c(0,54750))
par(new=T)
plot(t0_t12_rank_FPKM_Imp_lim_re$ourT12WTavg, col="green", pch =16, cex = 0.8, xlim=c(0,40), ylim=c(0,54750))
par(new=T)
plot(t0_t12_rank_FPKM_Imp_lim_re$ourT12ZFP57KOavg, col="red", pch =16, cex = 0.8, xlim=c(0,40), ylim=c(0,54750))


#Rank change
t0_t12_rank_FPKM_Imp_lim_re2 <- t0_t12_rank_FPKM_Imp_lim_re
t0_t12_rank_FPKM_Imp_lim_re2["ratioT0"] <- t0_t12_rank_FPKM_Imp_lim_re2$ourT0ZFP57KOJB1 - t0_t12_rank_FPKM_Imp_lim_re2$ourT0WTJB1

t0_t12_rank_FPKM_Imp_lim_re2["ratioT12"]  <- t0_t12_rank_FPKM_Imp_lim_re2$ourT12ZFP57KOavg -t0_t12_rank_FPKM_Imp_lim_re2$ourT12WTavg

plot(t0_t12_rank_FPKM_Imp_lim_re2$ratioT0, col="blue", pch =1, cex = 0.8, xlim=c(0,40), ylim=c(0,5))
par(new=T)
plot(t0_t12_rank_FPKM_Imp_lim_re2$ratioT12, col="brown", pch =1, cex = 0.8, xlim=c(0,40), ylim=c(0,5))
lines(x = c(0,40), y = c(1,1), col = "#464647", lty = 2)
#Save manually as ratiot0_t12_rank_FPKM_Imp_lim_re.svg
t0_t12_rank_FPKM_Imp_lim_re2st <- stack(as.matrix(t0_t12_rank_FPKM_Imp_lim_re2[,9:10]))
colnames(t0_t12_rank_FPKM_Imp_lim_re2st) <- c("Gene", "Group", "Ratio")
t0_t12_rank_FPKM_Imp_lim_re2st <- data.frame(t0_t12_rank_FPKM_Imp_lim_re2st)
library(ggplot2)
ggplot(t0_t12_rank_FPKM_Imp_lim_re2st) + coord_flip()+ ylim(c(-5000,5000))+
  geom_point(aes(x = Gene, y = Ratio, fill = Group, color = Group))+
  scale_color_manual(values = c("blue","red")) 

#t0_t12_rank_FPKM_Imp_lim_re2st.svg
breaksList12 = seq(-2000, 2000, by = 10)
pheatmap(t0_t12_rank_FPKM_Imp_lim_re2[,9:10],
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList12)),
         breaks = breaksList12,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=F, 
         cluster_rows=F,
         cellwidth = 20, 
         cellheight = 12)
#t0_t12_rank_FPKM_Imp_lim_re2.svg
dev.off()

head(t0_t12_rank_FPKM_Imp_lim_re2)
write.table(t0_t12_rank_FPKM_Imp_lim_re2, "t0_t12_rank_FPKM_Imp_lim_re2.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)

#DE imprinted genes




cp /media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/Bulk_mm10/featurecount/Igf2/imprinted_deregulated_genes_deseq2_0.05.sorted_JB1_Zfp57_overlapped_Kap1_mm10_chrMY_dedup.geneslist.txt imprinted_de_deseq2_0.05.txt  

fgrep -f imprinted_de_deseq2_0.05.txt  ourT0WTJB1_FPKM.txt -w | grep "mt-Th" -v | sort -k2,2 > ourT0WTJB1_FPKM_deImp_lim.txt
fgrep -f imprinted_de_deseq2_0.05.txt  ourT0ZFP57KOJB1_FPKM.txt -w | grep "mt-Th" -v | sort -k2,2 > ourT0ZFP57KOJB1_FPKM_deImp_lim.txt
fgrep -f imprinted_de_deseq2_0.05.txt  ourT12WT1_FPKM.txt -w | grep "mt-Th" -v | sort -k2,2 > ourT12WT1_FPKM_deImp_lim.txt
fgrep -f imprinted_de_deseq2_0.05.txt  ourT12WT2_FPKM.txt -w | grep "mt-Th" -v | sort -k2,2 > ourT12WT2_FPKM_deImp_lim.txt
fgrep -f imprinted_de_deseq2_0.05.txt  ourT12ZFP57KO1_FPKM.txt -w | grep "mt-Th" -v | sort -k2,2 > ourT12ZFP57KO1_FPKM_deImp_lim.txt
fgrep -f imprinted_de_deseq2_0.05.txt  ourT12ZFP57KO2_FPKM.txt -w | grep "mt-Th" -v | sort -k2,2 > ourT12ZFP57KO2_FPKM_deImp_lim.txt

fgrep -f imprinted_de_deseq2_0.05.txt  ourT12WTavg_FPKM.txt -w | grep "mt-Th" -v | sort -k2,2 > ourT12WTavg_FPKM_deImp_lim.txt
fgrep -f imprinted_de_deseq2_0.05.txt  ourT12ZFP57KOavg_FPKM.txt -w | grep "mt-Th" -v | sort -k2,2 > ourT12ZFP57KOavg_FPKM_deImp_lim.txt

paste *deImp_lim.txt  > t0_t12_rank_FPKM_deImp_lim.txt


t0_t12_rank_FPKM_deImp_lim <- read.table("t0_t12_rank_FPKM_deImp_lim.txt", header = F, stringsAsFactors = F)
rownames(t0_t12_rank_FPKM_deImp_lim) <- t0_t12_rank_FPKM_deImp_lim$V2
head(t0_t12_rank_FPKM_deImp_lim,2)
t0_t12_rank_FPKM_deImp_lim_re <- t0_t12_rank_FPKM_deImp_lim[,c(3,7,11,15,19,23,27,31)]
head(t0_t12_rank_FPKM_deImp_lim_re)
colnames(t0_t12_rank_FPKM_deImp_lim_re) <- c("ourT0WTJB1","ourT0ZFP57KOJB1","ourT12WT1","ourT12WT2","ourT12WTavg","ourT12ZFP57KO1","ourT12ZFP57KO2","ourT12ZFP57KOavg")
head(t0_t12_rank_FPKM_deImp_lim_re)
dim(t0_t12_rank_FPKM_deImp_lim_re)
write.table(t0_t12_rank_FPKM_deImp_lim_re, "t0_t12_rank_FPKM_deImp_lim_re.txt", sep="\t", quote = FALSE, append = FALSE, row.names = T)
with(t0_t12_rank_FPKM_deImp_lim_re, plot(ourT0WTJB1, ourT0ZFP57KOJB1, pch=21, cex = 0.2, main="Scatter Plot", xlim=c(0,54750), ylim=c(0,54750),lwd = 0.3,col="#525254", bg="grey",ylab="JB1 Zfp57-/- Rank", xlab="JB1 WT Rank", bty = 'n'))
with(subset(t0_t12_rank_FPKM_deImp_lim_re), points(ourT0WTJB1, ourT0ZFP57KOJB1, pch=21, cex = 0.8, xlim=c(0,54750), ylim=c(0,54750), lwd = 0.5,col="blue", bg="blue"),ylab="", xlab="", bty = 'n') #reddish
with(subset(t0_t12_rank_FPKM_deImp_lim_re), points(ourT12WTavg, ourT12ZFP57KOavg, pch=21, cex = 0.8, xlim=c(0,54750), ylim=c(0,54750), lwd = 0.5,col="red", bg="red"),ylab="", xlab="", bty = 'n') #reddish

lines(x = c(0,54750), y = c(0,54750), col = "#464647", lty = 2)
#Save manually as t0_t12_rank_FPKM_deImp_lim_re.svg


plot(t0_t12_rank_FPKM_deImp_lim_re$ourT0WTJB1, col="darkgreen", pch =16, cex = 0.8, xlim=c(0,40), ylim=c(0,54750))
par(new=T)
plot(t0_t12_rank_FPKM_deImp_lim_re$ourT0ZFP57KOJB1, col="darkred", pch =16, cex = 0.8, xlim=c(0,40), ylim=c(0,54750))
par(new=T)
plot(t0_t12_rank_FPKM_deImp_lim_re$ourT12WTavg, col="green", pch =16, cex = 0.8, xlim=c(0,40), ylim=c(0,54750))
par(new=T)
plot(t0_t12_rank_FPKM_deImp_lim_re$ourT12ZFP57KOavg, col="red", pch =16, cex = 0.8, xlim=c(0,40), ylim=c(0,54750))


#Rank change
t0_t12_rank_FPKM_deImp_lim_re2 <- t0_t12_rank_FPKM_deImp_lim_re
t0_t12_rank_FPKM_deImp_lim_re2["ratioT0KOvsWT"] <- t0_t12_rank_FPKM_deImp_lim_re2$ourT0ZFP57KOJB1 - t0_t12_rank_FPKM_deImp_lim_re2$ourT0WTJB1

t0_t12_rank_FPKM_deImp_lim_re2["ratioT12KOvsWT"]  <- t0_t12_rank_FPKM_deImp_lim_re2$ourT12ZFP57KOavg - t0_t12_rank_FPKM_deImp_lim_re2$ourT12WTavg
t0_t12_rank_FPKM_deImp_lim_re2["ratioT12WTvsT0WT"] <- t0_t12_rank_FPKM_deImp_lim_re2$ourT12WTavg - t0_t12_rank_FPKM_deImp_lim_re2$ourT0WTJB1

t0_t12_rank_FPKM_deImp_lim_re2["ratioT12KOvsT0KO"]  <- t0_t12_rank_FPKM_deImp_lim_re2$ourT12ZFP57KOavg - t0_t12_rank_FPKM_deImp_lim_re2$ourT0ZFP57KOJB1



plot(t0_t12_rank_FPKM_deImp_lim_re2$ratioT0, col="blue", pch =1, cex = 0.8, xlim=c(0,40), ylim=c(0,32))
par(new=T)
plot(t0_t12_rank_FPKM_deImp_lim_re2$ratioT12, col="brown", pch =1, cex = 0.8, xlim=c(0,40), ylim=c(0,32))
lines(x = c(0,40), y = c(1,1), col = "#464647", lty = 2)
#Save manually as ratiot0_t12_rank_FPKM_deImp_lim_re.svg

pheatmap(t0_t12_rank_FPKM_deImp_lim_re2[,9:12],
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList12)),
         breaks = breaksList12,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cellwidth = 20, 
         cellheight = 12)
#t0_t12_rank_FPKM_deImp_lim_re2.svg
#t0_t12_rank_FPKM_deImp_lim_re2.svg.jpeg
dev.off()



#####################################  END OF ANALYSIS #########################























#Filter samples and filter genes
#fgrep -f ./../../2019/time12/time0-tim12/removegenescellmarkers.txt bhybcountdatauq.cellmarkers1.norm.sep.txt  -v > bhybcountdatauq.cellmarkers1.norm_rmvfilt.txt
# Compute distances and hierarchical clustering
bhybcountdatauq.cellmarkers1.norm_rmvfilt <- read.table("bhybcountdatauq.cellmarkers1.norm_rmvfilt.txt", header = T, row.names = 1)
dim(bhybcountdatauq.cellmarkers1.norm_rmvfilt)
head(bhybcountdatauq.cellmarkers1.norm_rmvfilt)
#Remove samples E14 cell types
Tbhybcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp <- t(bhybcountdatauq.cellmarkers1.norm_rmvfilt[,c(3:12,22:25,28:31)])
bhybcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp_scale <- scale(Tbhybcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp)

ddbhybcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp_scale <- dist(as.matrix(bhybcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp_scale), method = "euclidean")
ddbhybcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp_scalehc <- hclust(ddbhybcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp_scale, method = "ward.D")
plot(ddbhybcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp_scalehc)
#save as ddbhybcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp_scalehc.svg


#z-score convert
z_Tbhybcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp= scale(t(bhybcountdatauq.cellmarkers1.norm_rmvfilt[,c(3:12,22:25,28:31)]), center = TRUE, scale = TRUE)
z_bhybcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp <- t(z_Tbhybcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp)

head(z_bhybcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp)
dim(z_bhybcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp)
#filtered samples
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_bhybcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp,
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12,
         filename = "pheatmap_z_bhybcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp.png")
dev.off()

# Compute distances and hierarchical clustering
bhybcountdatauq.cellmarkers1.norm_rmvfilt <- read.table("bhybcountdatauq.cellmarkers1.norm_rmvfilt.txt", header = T, row.names = 1)
dim(bhybcountdatauq.cellmarkers1.norm_rmvfilt)
head(bhybcountdatauq.cellmarkers1.norm_rmvfilt)
##Basilia Suggestion:Hybrid Type1 category
#Remove samples E14 cell types, ferguson E14 cell types
Tbhybcountdatauq.cellmarkers1.norm_rmvfilt_type1 <- t(bhybcountdatauq.cellmarkers1.norm_rmvfilt[,c(3:12,22:25)])
bhybcountdatauq.cellmarkers1.norm_rmvfilt_type1_scale <- scale(Tbhybcountdatauq.cellmarkers1.norm_rmvfilt_type1)

ddbhybcountdatauq.cellmarkers1.norm_rmvfilt_type1_scale <- dist(as.matrix(bhybcountdatauq.cellmarkers1.norm_rmvfilt_type1_scale), method = "euclidean")
ddbhybcountdatauq.cellmarkers1.norm_rmvfilt_type1_scalehc <- hclust(ddbhybcountdatauq.cellmarkers1.norm_rmvfilt_type1_scale, method = "ward.D")
plot(ddbhybcountdatauq.cellmarkers1.norm_rmvfilt_type1_scalehc)
#save as ddbhybcountdatauq.cellmarkers1.norm_rmvfilt_type1_scalehc.svg


#z-score convert
z_Tbhybcountdatauq.cellmarkers1.norm_rmvfilt_type1= scale(t(bhybcountdatauq.cellmarkers1.norm_rmvfilt[,c(3:12,22:25)]), center = TRUE, scale = TRUE)
z_bhybcountdatauq.cellmarkers1.norm_rmvfilt_type1 <- t(z_Tbhybcountdatauq.cellmarkers1.norm_rmvfilt_type1)

head(z_bhybcountdatauq.cellmarkers1.norm_rmvfilt_type1)
dim(z_bhybcountdatauq.cellmarkers1.norm_rmvfilt_type1)
#filtered samples
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_bhybcountdatauq.cellmarkers1.norm_rmvfilt_type1,
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12,
         filename = "pheatmap_z_bhybcountdatauq.cellmarkers1.norm_rmvfilt_type1.png")
dev.off()

##Basilia Suggestion: type2 category
#keep samples E14 cell types, ferguson E14 cell types and in vitro E14 bouschet
Tbhybcountdatauq.cellmarkers1.norm_rmvfilt_type2 <- t(bhybcountdatauq.cellmarkers1.norm_rmvfilt[,c(1:15,22:25,28:31)])
bhybcountdatauq.cellmarkers1.norm_rmvfilt_type2_scale <- scale(Tbhybcountdatauq.cellmarkers1.norm_rmvfilt_type2)

ddbhybcountdatauq.cellmarkers1.norm_rmvfilt_type2_scale <- dist(as.matrix(bhybcountdatauq.cellmarkers1.norm_rmvfilt_type2_scale), method = "euclidean")
ddbhybcountdatauq.cellmarkers1.norm_rmvfilt_type2_scalehc <- hclust(ddbhybcountdatauq.cellmarkers1.norm_rmvfilt_type2_scale, method = "ward.D")
plot(ddbhybcountdatauq.cellmarkers1.norm_rmvfilt_type2_scalehc)
#save as ddbhybcountdatauq.cellmarkers1.norm_rmvfilt_type2_scalehc.svg
library(ggdendro)
library(ggplot2)
fviz_dend(ddbhybcountdatauq.cellmarkers1.norm_rmvfilt_type2_scalehc, 
          cex = 1, 
          k = 3, 
          k_colors = c("#D95F02", "#D95F02", "#6B8E23"),
          palette= c("#D95F02", "#D95F02", "#6B8E23"),
          horiz = TRUE,
          color_labels_by_k = TRUE,
          type = "rectangle",
          labels_track_height= 25)
#save as fvizddbhybcountdatauq.cellmarkers1.norm_rmvfilt_type2_scalehc.png

#z-score convert
z_Tbhybcountdatauq.cellmarkers1.norm_rmvfilt_type2= scale(t(bhybcountdatauq.cellmarkers1.norm_rmvfilt[,c(1:15,22:25,28:31)]), center = TRUE, scale = TRUE)
z_bhybcountdatauq.cellmarkers1.norm_rmvfilt_type2 <- t(z_Tbhybcountdatauq.cellmarkers1.norm_rmvfilt_type2)

head(z_bhybcountdatauq.cellmarkers1.norm_rmvfilt_type2)
dim(z_bhybcountdatauq.cellmarkers1.norm_rmvfilt_type2)
#filtered samples
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_bhybcountdatauq.cellmarkers1.norm_rmvfilt_type2,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(length(breaksList1)),
         breaks = breaksList1,
         treeheight_row = 30,
         treeheight_col = 30,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 10, 
         cellheight = 10, filename = "pheatmap_z_bhybcountdatauq.cellmarkers1.norm_rmvfilt_type2.png")
dev.off()



###########################################################################à

#Bulk
#folder in cluster 
#OurT0
/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotationlike-_chr.gtf -T 12 -o Bulk_ourT0_star-featureCounts_GRCm38.mm10.txt ourT0_WT_Rep1_Aligned.sortedByCoord.out.bam ourT0_WT_Rep2_Aligned.sortedByCoord.out.bam ourT0_WT_Rep3_Aligned.sortedByCoord.out.bam ourT0_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bam ourT0_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bam ourT0_ZFP57_KO_Rep3_Aligned.sortedByCoord.out.bam
#Our T0 Riso GSE77444, E14 n JB1
/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotationlike-_chr.gtf -T 12 -o Bulk_ourT0ej_star-featureCounts_GRCm38.mm10.txt SRR3137616_Aligned.sortedByCoord.out.bam SRR3137617_Aligned.sortedByCoord.out.bam SRR3137618_Aligned.sortedByCoord.out.bam SRR3137619_Aligned.sortedByCoord.out.bam
#Bouschet
/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotationlike-_chr.gtf -T 12 -o Bulk_bouchet_star-featureCounts_GRCm38.mm10.txt SRR1409932_Aligned.sortedByCoord.out.bam SRR1409933_Aligned.sortedByCoord.out.bam SRR1409934_Aligned.sortedByCoord.out.bam SRR1409935_Aligned.sortedByCoord.out.bam SRR1409938_Aligned.sortedByCoord.out.bam SRR1409939_Aligned.sortedByCoord.out.bam SRR1409940_Aligned.sortedByCoord.out.bam SRR1409941_Aligned.sortedByCoord.out.bam SRR1409944_Aligned.sortedByCoord.out.bam SRR1409945_Aligned.sortedByCoord.out.bam SRR1409946_Aligned.sortedByCoord.out.bam SRR1409950_Aligned.sortedByCoord.out.bam SRR1409951_Aligned.sortedByCoord.out.bam SRR1409952_Aligned.sortedByCoord.out.bam SRR1409953_Aligned.sortedByCoord.out.bam SRR1409954_Aligned.sortedByCoord.out.bam SRR1409955_Aligned.sortedByCoord.out.bam
#T12
/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -p -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotationlike-_chr.gtf -T 12 -o Bulk_ourT12_star-featureCounts_GRCm38.mm10.txt ourT12_WT_Rep1_Aligned.sortedByCoord.out.bam ourT12_WT_Rep2_Aligned.sortedByCoord.out.bam ourT12_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bam ourT12_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bam
#Feil
/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotationlike-_chr.gtf -T 12 -o Bulk_wtNPC_feil_star-featureCounts_mm10_Igf2.txt SRR5665889_Aligned.sortedByCoord.out.bam SRR5665890_Aligned.sortedByCoord.out.bam
#Ferguson
/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -p -s 1 -a /home/ankitv/ref_av/gencodes/gencode_M20/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotationlike-_chr.gtf -T 12 -o Bulk_esc_ferg_star-featureCounts_GRCm38.mm10.txt SRR8329324_Aligned.sortedByCoord.out.bam SRR8329325_Aligned.sortedByCoord.out.bam SRR8329327_Aligned.sortedByCoord.out.bam SRR8329328_Aligned.sortedByCoord.out.bam


################################################################# DESeq2 ######################################################################
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2020/CombineOurNBouschetNFeilNFerguson/")
bcountdata1 <- read.table("Bulk_ourT0ej_star-featureCounts_GRCm38.mm10.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
bcountdata1 <- bcountdata1[ ,6:ncol(bcountdata1)]
head(bcountdata1)
# Remove .bam or .sam from filenames
colnames(bcountdata1) <- gsub("\\.[sb]am$", "", colnames(bcountdata1))
colnames(bcountdata1) <-  c("ourT0WTE14","ourT0ZFP57KOE14","ourT0WTJB1","ourT0ZFP57KOJB1")
bcountdata1 <- bcountdata1[order(rownames(bcountdata1)),]
head(bcountdata1)
dim(bcountdata1)
tail(bcountdata1)
library(NOISeq)
bcountdata1uq <- uqua(bcountdata1, long = 1000, lc = 0, k = 0)
head(bcountdata1uq)


bcountdata2 <- read.table("Bulk_bouchet_star-featureCounts_GRCm38.mm10.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
bcountdata2 <- bcountdata2[ ,6:ncol(bcountdata2)]
head(bcountdata2)
# Remove .bam or .sam from filenames
colnames(bcountdata2) <- gsub("\\.[sb]am$", "", colnames(bcountdata2))
colnames(bcountdata2) <-  c("SRR1409932","SRR1409933","SRR1409934","SRR1409935","SRR1409938","SRR1409939","SRR1409940","SRR1409941","SRR1409944","SRR1409945","SRR1409946","SRR1409950","SRR1409951","SRR1409952","SRR1409953","SRR1409954","SRR1409955")
bcountdata2 <- bcountdata2[order(rownames(bcountdata2)),]
head(bcountdata2)
dim(bcountdata2)
tail(bcountdata2)
bcountdata2uq <- uqua(bcountdata2, long = 1000, lc = 0, k = 0)
head(bcountdata2uq)

bcountdata3 <- read.table("Bulk_ourT12_star-featureCounts_GRCm38.mm10.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
bcountdata3 <- bcountdata3[ ,6:ncol(bcountdata3)]
head(bcountdata3)
# Remove .bam or .sam from filenames
colnames(bcountdata3) <- gsub("\\.[sb]am$", "", colnames(bcountdata3))
colnames(bcountdata3) <-  c("ourT12WT1","ourT12WT2","ourT12ZFP57KO1","ourT12ZFP57KO2")
bcountdata3 <- bcountdata3[order(rownames(bcountdata3)),]
head(bcountdata3)
dim(bcountdata3)
tail(bcountdata3)
bcountdata3uq <- uqua(bcountdata3, long = 1000, lc = 0, k = 0)
head(bcountdata3uq)

bcountdata4 <- read.table("Bulk_wtNPC_feil_star-featureCounts_mm10_Igf2.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
bcountdata4 <- bcountdata4[ ,6:ncol(bcountdata4)]
head(bcountdata4)
# Remove .bam or .sam from filenames
colnames(bcountdata4) <- gsub("\\.[sb]am$", "", colnames(bcountdata4))
colnames(bcountdata4) <-  c("SRR5665889","SRR5665890")
bcountdata4 <- bcountdata4[order(rownames(bcountdata4)),]
head(bcountdata4)
dim(bcountdata4)
tail(bcountdata4)
bcountdata4uq <- uqua(bcountdata4, long = 1000, lc = 0, k = 0)
head(bcountdata4uq)

bcountdata5 <- read.table("Bulk_esc_ferg_star-featureCounts_GRCm38.mm10.txt", header=TRUE, row.names=1)
# Remove first five columns (chr, start, end, strand, length)
bcountdata5 <- bcountdata5[ ,6:ncol(bcountdata5)]
head(bcountdata5)
# Remove .bam or .sam from filenames
colnames(bcountdata5) <- gsub("\\.[sb]am$", "", colnames(bcountdata5))
colnames(bcountdata5) <-  c("SRR8329324","SRR8329325","SRR8329327","SRR8329328")
bcountdata5 <- bcountdata5[order(rownames(bcountdata5)),]
head(bcountdata5)
dim(bcountdata5)
tail(bcountdata5)
bcountdata5uq <- uqua(bcountdata5, long = 1000, lc = 0, k = 0)
head(bcountdata5uq)

bcountdata <- cbind.data.frame(bcountdata1, bcountdata2, bcountdata3, bcountdata4, bcountdata5)
library(DESeq2)
bcountdata = as.matrix(bcountdata)
head(bcountdata)
dim(bcountdata)
boxplot(bcountdata, ylim=c(0,100))



bcountdatauq <- cbind.data.frame(bcountdata1uq, bcountdata2uq, bcountdata3uq, bcountdata4uq, bcountdata5uq)
library(DESeq2)
bcountdatauq = data.frame(bcountdatauq)
head(bcountdatauq)
dim(bcountdatauq)

#boxplot(bcountdatauq, ylim=c(0,100))
bcountdatauq["id"] <- rownames(bcountdatauq)

#sort -k1,1 -u cellmarkerslist.txt | sort -k2,2 > cellmarkerslist.sorted.txt
cellmarkerslist.sorted <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/cellmarkerslist.sorted.txt", header = F)
head(cellmarkerslist.sorted)
colnames(cellmarkerslist.sorted) <- c("Gene", "markers")
dim(cellmarkerslist.sorted) #print 94 genes
ens_gene_names_chrpos_dedup_M20 <- read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt")
head(ens_gene_names_chrpos_dedup_M20)
ens_gene_names_chrpos_dedup_M20 <- ens_gene_names_chrpos_dedup_M20[,c(3,4,5,1,2)]
head(ens_gene_names_chrpos_dedup_M20)
colnames(ens_gene_names_chrpos_dedup_M20) <- c("chr", "start", "end", "id", "Gene")

cellmarkers = merge(cellmarkerslist.sorted, ens_gene_names_chrpos_dedup_M20, by="Gene", all.x=FALSE)
head(cellmarkers)
dim(cellmarkers) #print all 94 genes assigned
cellmarkersre <- cellmarkers[order(cellmarkers$markers),]
head(cellmarkersre)
dim(cellmarkersre)
#cellmarkers.cords <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/cellmarkers.symbol.ens.chr.txt", header = F)
#head(cellmarkers.cords)
#dim(cellmarkers.cords)
#colnames(cellmarkers.cords) <- c("chr", "start", "end", "id", "Gene")

bcountdatauq.cellmarker = merge(bcountdatauq, cellmarkersre, by="id", all.x=FALSE)
head(bcountdatauq.cellmarker)
dim(bcountdatauq.cellmarker) #Print all 94 genes overlapped and assigned
bcountdatauq.cellmarkers <- bcountdatauq.cellmarker[order(bcountdatauq.cellmarker$markers),]
dim(bcountdatauq.cellmarkers)
head(bcountdatauq.cellmarkers)
bcountdatauq.cellmarkersre <- bcountdatauq.cellmarkers[,c(33, 2:32)]
head(bcountdatauq.cellmarkersre,2)
rownames(bcountdatauq.cellmarkersre) <- bcountdatauq.cellmarkersre[,1]
bcountdatauq.cellmarkers1 <- bcountdatauq.cellmarkersre[,c(2:32)]
head(bcountdatauq.cellmarkers1,2)
bcountdatauq.cellmarkers1 <- as.matrix(bcountdatauq.cellmarkers1)

#Loading control
loadingcontrol.cords <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/loadingcontrol.symbol.ens.chr.txt", header = F)
head(loadingcontrol.cords)
colnames(loadingcontrol.cords) <- c("chr", "start", "end", "id", "Gene")
bcountdatauq.loadingcontrol = merge(bcountdatauq, loadingcontrol.cords, by="id", all.x=FALSE)
head(bcountdatauq.loadingcontrol)
dim(bcountdatauq.loadingcontrol)
bcountdatauq.loadingcontrol <- bcountdatauq.loadingcontrol[,c(33, 2:32)]
head(bcountdatauq.loadingcontrol,2)
rownames(bcountdatauq.loadingcontrol) <- bcountdatauq.loadingcontrol$Gene
bcountdatauq.loadingcontrol1 <- bcountdatauq.loadingcontrol[,c(2:32)]
head(bcountdatauq.loadingcontrol1,2)
bcountdatauq.loadingcontrol1 <- as.matrix(bcountdatauq.loadingcontrol1)
bcountdatauq.loadingcontrol1.Actb <- bcountdatauq.loadingcontrol1[1,]
#Normalize with loading control
dim(bcountdatauq.cellmarkers1)

bcountdatauq.cellmarkers1.norm <- cbind.data.frame((bcountdatauq.cellmarkers1[,1] / bcountdatauq.loadingcontrol1.Actb[1]),
                                                   (bcountdatauq.cellmarkers1[,2] / bcountdatauq.loadingcontrol1.Actb[2]),
                                                   (bcountdatauq.cellmarkers1[,3] / bcountdatauq.loadingcontrol1.Actb[3]),
                                                   (bcountdatauq.cellmarkers1[,4] / bcountdatauq.loadingcontrol1.Actb[4]),
                                                   (bcountdatauq.cellmarkers1[,5] / bcountdatauq.loadingcontrol1.Actb[5]),
                                                   (bcountdatauq.cellmarkers1[,6] / bcountdatauq.loadingcontrol1.Actb[6]),
                                                   (bcountdatauq.cellmarkers1[,7] / bcountdatauq.loadingcontrol1.Actb[7]),
                                                   (bcountdatauq.cellmarkers1[,8] / bcountdatauq.loadingcontrol1.Actb[8]),
                                                   (bcountdatauq.cellmarkers1[,9] / bcountdatauq.loadingcontrol1.Actb[9]),
                                                   (bcountdatauq.cellmarkers1[,10] / bcountdatauq.loadingcontrol1.Actb[10]),
                                                   (bcountdatauq.cellmarkers1[,11] / bcountdatauq.loadingcontrol1.Actb[11]),
                                                   (bcountdatauq.cellmarkers1[,12] / bcountdatauq.loadingcontrol1.Actb[12]),
                                                   (bcountdatauq.cellmarkers1[,13] / bcountdatauq.loadingcontrol1.Actb[13]),
                                                   (bcountdatauq.cellmarkers1[,14] / bcountdatauq.loadingcontrol1.Actb[14]),
                                                   (bcountdatauq.cellmarkers1[,15] / bcountdatauq.loadingcontrol1.Actb[15]),
                                                   (bcountdatauq.cellmarkers1[,16] / bcountdatauq.loadingcontrol1.Actb[16]),
                                                   (bcountdatauq.cellmarkers1[,17] / bcountdatauq.loadingcontrol1.Actb[17]),
                                                   (bcountdatauq.cellmarkers1[,18] / bcountdatauq.loadingcontrol1.Actb[18]),
                                                   (bcountdatauq.cellmarkers1[,19] / bcountdatauq.loadingcontrol1.Actb[19]),
                                                   (bcountdatauq.cellmarkers1[,20] / bcountdatauq.loadingcontrol1.Actb[20]),
                                                   (bcountdatauq.cellmarkers1[,21] / bcountdatauq.loadingcontrol1.Actb[21]),
                                                   (bcountdatauq.cellmarkers1[,22] / bcountdatauq.loadingcontrol1.Actb[22]),
                                                   (bcountdatauq.cellmarkers1[,23] / bcountdatauq.loadingcontrol1.Actb[23]),
                                                   (bcountdatauq.cellmarkers1[,24] / bcountdatauq.loadingcontrol1.Actb[24]),
                                                   (bcountdatauq.cellmarkers1[,25] / bcountdatauq.loadingcontrol1.Actb[25]),
                                                   (bcountdatauq.cellmarkers1[,26] / bcountdatauq.loadingcontrol1.Actb[26]),
                                                   (bcountdatauq.cellmarkers1[,27] / bcountdatauq.loadingcontrol1.Actb[27]),
                                                   (bcountdatauq.cellmarkers1[,28] / bcountdatauq.loadingcontrol1.Actb[28]),
                                                   (bcountdatauq.cellmarkers1[,29] / bcountdatauq.loadingcontrol1.Actb[29]),
                                                   (bcountdatauq.cellmarkers1[,30] / bcountdatauq.loadingcontrol1.Actb[30]),
                                                   (bcountdatauq.cellmarkers1[,31] / bcountdatauq.loadingcontrol1.Actb[31]))
colnames(bcountdatauq.cellmarkers1.norm) <- c("ourT0WTE14","ourT0ZFP57KOE14","ourT0WTJB1","ourT0ZFP57KOJB1","SRR1409932","SRR1409933","SRR1409934","SRR1409935","SRR1409938","SRR1409939","SRR1409940","SRR1409941","SRR1409944","SRR1409945","SRR1409946","SRR1409950","SRR1409951","SRR1409952","SRR1409953","SRR1409954","SRR1409955","ourT12WT1", "ourT12WT2", "ourT12ZFP57KO1", "ourT12ZFP57KO2","SRR5665889","SRR5665890","SRR8329324","SRR8329325","SRR8329327","SRR8329328")

head(bcountdatauq.cellmarkers1.norm)
boxplot(log10(bcountdatauq.cellmarkers1.norm), ylim=c(0,10), col = c("blue","blue","blue","blue","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","grey","grey","orange","orange","orange","orange"))
#save as bcountdatauq.cellmarkers1.norm.sep.svg
library(EDASeq)
plotPCA(as.matrix(bcountdatauq.cellmarkers1.norm))
plotPCA(as.matrix(bcountdatauq.cellmarkers1.norm), labels=F, col =  c("cyan","red","cyan","red","blue","blue","darkgreen","darkgreen","blue","blue","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","seagreen","seagreen","seagreen","seagreen","seagreen","seagreen","green","green","magenta","magenta","darkgreen","darkgreen","navy","navy","red","red"))
ggsave("PCA_tbcountdatauq.cellmarkers1.norm_scaleT.sep.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

write.table(bcountdatauq.cellmarkers1.norm, "bcountdatauq.cellmarkers1.norm.sep.txt", row.names = T, quote = F, append=F)
#Dendrogram cluster
dim(bcountdatauq.cellmarkers1.norm)  #94 31
# Compute distances and hierarchical clustering
bcountdatauq.cellmarkers1.norm1 <- t(bcountdatauq.cellmarkers1.norm)
bcountdatauq.cellmarkers1.normscale <- scale(bcountdatauq.cellmarkers1.norm1)
ddbcountdatauq.cellmarkers1.normscale <- dist(as.matrix(bcountdatauq.cellmarkers1.normscale), method = "euclidean")
ddbcountdatauq.cellmarkers1.normscalehc <- hclust(ddbcountdatauq.cellmarkers1.normscale, method = "ward.D")
plot(ddbcountdatauq.cellmarkers1.normscalehc)
#save as ddbcountdatauq.cellmarkers1.normscalehc.sep.svg
head(bcountdatauq.cellmarkers1.norm,1)                                                             
colnames(bcountdatauq.cellmarkers1.norm) <- colnames(bcountdatauq.cellmarkers1)
head(bcountdatauq.cellmarkers1.norm, 2)                                                             
z_Tbcountdatauq.cellmarkers1.norm= scale(t(bcountdatauq.cellmarkers1.norm), center = TRUE, scale = TRUE)
z_bcountdatauq.cellmarkers1.norm <- t(z_Tbcountdatauq.cellmarkers1.norm)
colfunc <- colorRampPalette(c("navy","white","firebrick3"))

dim(z_bcountdatauq.cellmarkers1.norm)  #94 31
library(pheatmap)
library(RColorBrewer)
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_bcountdatauq.cellmarkers1.norm[,c(1:31)],
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12, 
         filename = "pheatmap_z_bcountdatauq.cellmarkers1.norm.sep.png")

dev.off()
#Dendrogram cluster
dim(bcountdatauq.cellmarkers1.norm)
head(bcountdatauq.cellmarkers1.norm)
#Filter samples and filter genes
#fgrep -f ./../../2019/time12/time0-tim12/removegenescellmarkers.txt bcountdatauq.cellmarkers1.norm.sep.txt  -v > bcountdatauq.cellmarkers1.norm_rmvfilt.txt
# Compute distances and hierarchical clustering
bcountdatauq.cellmarkers1.norm_rmvfilt = bcountdatauq.cellmarkers1.norm
bcountdatauq.cellmarkers1.norm_rmvfilt <- read.table("bcountdatauq.cellmarkers1.norm_rmvfilt.txt", header = T, row.names = 1)
dim(bcountdatauq.cellmarkers1.norm_rmvfilt)
head(bcountdatauq.cellmarkers1.norm_rmvfilt)
#Remove samples E14 cell types
Tbcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp <- t(bcountdatauq.cellmarkers1.norm_rmvfilt[,c(3:12,22:25,28:31)])
bcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp_scale <- scale(Tbcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp)

ddbcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp_scale <- dist(as.matrix(bcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp_scale), method = "euclidean")
ddbcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp_scalehc <- hclust(ddbcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp_scale, method = "ward.D")
plot(ddbcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp_scalehc)
#save as ddbcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp_scalehc.svg


#z-score convert
z_Tbcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp= scale(t(bcountdatauq.cellmarkers1.norm_rmvfilt[,c(3:12,22:25,28:31)]), center = TRUE, scale = TRUE)
z_bcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp <- t(z_Tbcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp)

head(z_bcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp)
dim(z_bcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp)
#filtered samples
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_bcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp,
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12,
         filename = "pheatmap_z_bcountdatauq.cellmarkers1.norm_rmvfilt_rmvsamp.png")
dev.off()
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#created after basila gene selection  bouschetmarkerslist.sorted.txt
bouschetmarkerslist.sorted <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/bouschetmarkerslist.sorted.txt", header = F)
head(bouschetmarkerslist.sorted)
colnames(bouschetmarkerslist.sorted) <- c("Gene", "markers")
dim(bouschetmarkerslist.sorted) #print 29 genes
ens_gene_names_chrpos_dedup_M20 <- read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt")
head(ens_gene_names_chrpos_dedup_M20)
ens_gene_names_chrpos_dedup_M20 <- ens_gene_names_chrpos_dedup_M20[,c(3,4,5,1,2)]
head(ens_gene_names_chrpos_dedup_M20)
colnames(ens_gene_names_chrpos_dedup_M20) <- c("chr", "start", "end", "id", "Gene")

bouschetmarkers = merge(bouschetmarkerslist.sorted, ens_gene_names_chrpos_dedup_M20, by="Gene", all.x=FALSE)
head(bouschetmarkers)
dim(bouschetmarkers) #print all 29 genes assigned
bouschetmarkersre <- bouschetmarkers[order(bouschetmarkers$markers),]
head(bouschetmarkersre)
dim(bouschetmarkersre)
#bouschetmarkers.cords <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/bouschetmarkers.symbol.ens.chr.txt", header = F)
#head(bouschetmarkers.cords)
#dim(bouschetmarkers.cords)
#colnames(bouschetmarkers.cords) <- c("chr", "start", "end", "id", "Gene")

bcountdatauq.bouschetmarker = merge(bcountdatauq, bouschetmarkersre, by="id", all.x=FALSE)
head(bcountdatauq.bouschetmarker)
dim(bcountdatauq.bouschetmarker) #Print all 29 genes overlapped and assigned
bcountdatauq.bouschetmarkers <- bcountdatauq.bouschetmarker[order(bcountdatauq.bouschetmarker$markers),]
dim(bcountdatauq.bouschetmarkers)
head(bcountdatauq.bouschetmarkers)
bcountdatauq.bouschetmarkersre <- bcountdatauq.bouschetmarkers[,c(33, 2:32)]
head(bcountdatauq.bouschetmarkersre,2)
rownames(bcountdatauq.bouschetmarkersre) <- bcountdatauq.bouschetmarkersre[,1]
bcountdatauq.bouschetmarkers1 <- bcountdatauq.bouschetmarkersre[,c(2:32)]
head(bcountdatauq.bouschetmarkers1,2)
bcountdatauq.bouschetmarkers1 <- as.matrix(bcountdatauq.bouschetmarkers1)

#Loading control
loadingcontrol.cords <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/loadingcontrol.symbol.ens.chr.txt", header = F)
head(loadingcontrol.cords)
colnames(loadingcontrol.cords) <- c("chr", "start", "end", "id", "Gene")
bcountdatauq.loadingcontrol = merge(bcountdatauq, loadingcontrol.cords, by="id", all.x=FALSE)
head(bcountdatauq.loadingcontrol)
dim(bcountdatauq.loadingcontrol)
bcountdatauq.loadingcontrol <- bcountdatauq.loadingcontrol[,c(33, 2:32)]
head(bcountdatauq.loadingcontrol,2)
rownames(bcountdatauq.loadingcontrol) <- bcountdatauq.loadingcontrol$Gene
bcountdatauq.loadingcontrol1 <- bcountdatauq.loadingcontrol[,c(2:32)]
head(bcountdatauq.loadingcontrol1,2)
bcountdatauq.loadingcontrol1 <- as.matrix(bcountdatauq.loadingcontrol1)
bcountdatauq.loadingcontrol1.Actb <- bcountdatauq.loadingcontrol1[1,]
#Normalize with loading control
dim(bcountdatauq.bouschetmarkers1)

bcountdatauq.bouschetmarkers1.norm <- cbind.data.frame((bcountdatauq.bouschetmarkers1[,1] / bcountdatauq.loadingcontrol1.Actb[1]),
                                                       (bcountdatauq.bouschetmarkers1[,2] / bcountdatauq.loadingcontrol1.Actb[2]),
                                                       (bcountdatauq.bouschetmarkers1[,3] / bcountdatauq.loadingcontrol1.Actb[3]),
                                                       (bcountdatauq.bouschetmarkers1[,4] / bcountdatauq.loadingcontrol1.Actb[4]),
                                                       (bcountdatauq.bouschetmarkers1[,5] / bcountdatauq.loadingcontrol1.Actb[5]),
                                                       (bcountdatauq.bouschetmarkers1[,6] / bcountdatauq.loadingcontrol1.Actb[6]),
                                                       (bcountdatauq.bouschetmarkers1[,7] / bcountdatauq.loadingcontrol1.Actb[7]),
                                                       (bcountdatauq.bouschetmarkers1[,8] / bcountdatauq.loadingcontrol1.Actb[8]),
                                                       (bcountdatauq.bouschetmarkers1[,9] / bcountdatauq.loadingcontrol1.Actb[9]),
                                                       (bcountdatauq.bouschetmarkers1[,10] / bcountdatauq.loadingcontrol1.Actb[10]),
                                                       (bcountdatauq.bouschetmarkers1[,11] / bcountdatauq.loadingcontrol1.Actb[11]),
                                                       (bcountdatauq.bouschetmarkers1[,12] / bcountdatauq.loadingcontrol1.Actb[12]),
                                                       (bcountdatauq.bouschetmarkers1[,13] / bcountdatauq.loadingcontrol1.Actb[13]),
                                                       (bcountdatauq.bouschetmarkers1[,14] / bcountdatauq.loadingcontrol1.Actb[14]),
                                                       (bcountdatauq.bouschetmarkers1[,15] / bcountdatauq.loadingcontrol1.Actb[15]),
                                                       (bcountdatauq.bouschetmarkers1[,16] / bcountdatauq.loadingcontrol1.Actb[16]),
                                                       (bcountdatauq.bouschetmarkers1[,17] / bcountdatauq.loadingcontrol1.Actb[17]),
                                                       (bcountdatauq.bouschetmarkers1[,18] / bcountdatauq.loadingcontrol1.Actb[18]),
                                                       (bcountdatauq.bouschetmarkers1[,19] / bcountdatauq.loadingcontrol1.Actb[19]),
                                                       (bcountdatauq.bouschetmarkers1[,20] / bcountdatauq.loadingcontrol1.Actb[20]),
                                                       (bcountdatauq.bouschetmarkers1[,21] / bcountdatauq.loadingcontrol1.Actb[21]),
                                                       (bcountdatauq.bouschetmarkers1[,22] / bcountdatauq.loadingcontrol1.Actb[22]),
                                                       (bcountdatauq.bouschetmarkers1[,23] / bcountdatauq.loadingcontrol1.Actb[23]),
                                                       (bcountdatauq.bouschetmarkers1[,24] / bcountdatauq.loadingcontrol1.Actb[24]),
                                                       (bcountdatauq.bouschetmarkers1[,25] / bcountdatauq.loadingcontrol1.Actb[25]),
                                                       (bcountdatauq.bouschetmarkers1[,26] / bcountdatauq.loadingcontrol1.Actb[26]),
                                                       (bcountdatauq.bouschetmarkers1[,27] / bcountdatauq.loadingcontrol1.Actb[27]),
                                                       (bcountdatauq.bouschetmarkers1[,28] / bcountdatauq.loadingcontrol1.Actb[28]),
                                                       (bcountdatauq.bouschetmarkers1[,29] / bcountdatauq.loadingcontrol1.Actb[29]),
                                                       (bcountdatauq.bouschetmarkers1[,30] / bcountdatauq.loadingcontrol1.Actb[30]),
                                                       (bcountdatauq.bouschetmarkers1[,31] / bcountdatauq.loadingcontrol1.Actb[31]))
colnames(bcountdatauq.bouschetmarkers1.norm) <- c("ourT0WTE14","ourT0ZFP57KOE14","ourT0WTJB1","ourT0ZFP57KOJB1","SRR1409932","SRR1409933","SRR1409934","SRR1409935","SRR1409938","SRR1409939","SRR1409940","SRR1409941","SRR1409944","SRR1409945","SRR1409946","SRR1409950","SRR1409951","SRR1409952","SRR1409953","SRR1409954","SRR1409955","ourT12WT1", "ourT12WT2", "ourT12ZFP57KO1", "ourT12ZFP57KO2","SRR5665889","SRR5665890","SRR8329324","SRR8329325","SRR8329327","SRR8329328")

head(bcountdatauq.bouschetmarkers1.norm)
boxplot(log10(bcountdatauq.bouschetmarkers1.norm), ylim=c(0,10), col = c("blue","blue","blue","blue","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","grey","grey","orange","orange","orange","orange"))
#save as bcountdatauq.bouschetmarkers1.norm.sep.svg
library(EDASeq)
plotPCA(as.matrix(bcountdatauq.bouschetmarkers1.norm))
plotPCA(as.matrix(bcountdatauq.bouschetmarkers1.norm), labels=F, col =  c("cyan","red","cyan","red","blue","blue","darkgreen","darkgreen","blue","blue","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","seagreen","seagreen","seagreen","seagreen","seagreen","seagreen","green","green","magenta","magenta","darkgreen","darkgreen","navy","navy","red","red"))
ggsave("PCA_tbcountdatauq.bouschetmarkers1.norm_scaleT.sep.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

write.table(bcountdatauq.bouschetmarkers1.norm, "bcountdatauq.bouschetmarkers1.norm.sep.txt", row.names = T, quote = F, append=F)
#Dendrogram cluster
dim(bcountdatauq.bouschetmarkers1.norm)  #29 31
# Compute distances and hierarchical clustering
bcountdatauq.bouschetmarkers1.norm1 <- t(bcountdatauq.bouschetmarkers1.norm)
bcountdatauq.bouschetmarkers1.normscale <- scale(bcountdatauq.bouschetmarkers1.norm1)
ddbcountdatauq.bouschetmarkers1.normscale <- dist(as.matrix(bcountdatauq.bouschetmarkers1.normscale), method = "euclidean")
ddbcountdatauq.bouschetmarkers1.normscalehc <- hclust(ddbcountdatauq.bouschetmarkers1.normscale, method = "ward.D")
plot(ddbcountdatauq.bouschetmarkers1.normscalehc)
#save as ddbcountdatauq.bouschetmarkers1.normscalehc.sep.svg
head(bcountdatauq.bouschetmarkers1.norm,1)                                                             
colnames(bcountdatauq.bouschetmarkers1.norm) <- colnames(bcountdatauq.bouschetmarkers1)
head(bcountdatauq.bouschetmarkers1.norm, 2)                                                             
z_Tbcountdatauq.bouschetmarkers1.norm= scale(t(bcountdatauq.bouschetmarkers1.norm), center = TRUE, scale = TRUE)
z_bcountdatauq.bouschetmarkers1.norm <- t(z_Tbcountdatauq.bouschetmarkers1.norm)
colfunc <- colorRampPalette(c("navy","white","firebrick3"))

dim(z_bcountdatauq.bouschetmarkers1.norm)  #29 31
library(pheatmap)
library(RColorBrewer)
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_bcountdatauq.bouschetmarkers1.norm[,c(1:31)],
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12, 
         filename = "pheatmap_z_bcountdatauq.bouschetmarkers1.norm.sep.png")

dev.off()
#Dendrogram cluster
dim(bcountdatauq.bouschetmarkers1.norm)
head(bcountdatauq.bouschetmarkers1.norm)


#Filter samples and filter genes
#fgrep -f ./../../2019/time12/time0-tim12/removegenesbouschetmarkers.txt bcountdatauq.bouschetmarkers1.norm.sep.txt  -v > bcountdatauq.bouschetmarkers1.norm_rmvfilt.txt
# Compute distances and hierarchical clustering
bcountdatauq.bouschetmarkers1.norm_rmvfilt <- read.table("bcountdatauq.bouschetmarkers1.norm_rmvfilt.txt", header = T, row.names = 1)
dim(bcountdatauq.bouschetmarkers1.norm_rmvfilt)
head(bcountdatauq.bouschetmarkers1.norm_rmvfilt)
#Remove samples E14 cell types
Tbcountdatauq.bouschetmarkers1.norm_rmvfilt_rmvsamp <- t(bcountdatauq.bouschetmarkers1.norm_rmvfilt[,c(3:12,22:25,28:31)])
bcountdatauq.bouschetmarkers1.norm_rmvfilt_rmvsamp_scale <- scale(Tbcountdatauq.bouschetmarkers1.norm_rmvfilt_rmvsamp)

ddbcountdatauq.bouschetmarkers1.norm_rmvfilt_rmvsamp_scale <- dist(as.matrix(bcountdatauq.bouschetmarkers1.norm_rmvfilt_rmvsamp_scale), method = "euclidean")
ddbcountdatauq.bouschetmarkers1.norm_rmvfilt_rmvsamp_scalehc <- hclust(ddbcountdatauq.bouschetmarkers1.norm_rmvfilt_rmvsamp_scale, method = "ward.D")
plot(ddbcountdatauq.bouschetmarkers1.norm_rmvfilt_rmvsamp_scalehc)
#save as ddbcountdatauq.bouschetmarkers1.norm_rmvfilt_rmvsamp_scalehc.svg


#z-score convert
z_Tbcountdatauq.bouschetmarkers1.norm_rmvfilt_rmvsamp= scale(t(bcountdatauq.bouschetmarkers1.norm_rmvfilt[,c(3:12,22:25,28:31)]), center = TRUE, scale = TRUE)
z_bcountdatauq.bouschetmarkers1.norm_rmvfilt_rmvsamp <- t(z_Tbcountdatauq.bouschetmarkers1.norm_rmvfilt_rmvsamp)

head(z_bcountdatauq.bouschetmarkers1.norm_rmvfilt_rmvsamp)
dim(z_bcountdatauq.bouschetmarkers1.norm_rmvfilt_rmvsamp)
#filtered samples
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_bcountdatauq.bouschetmarkers1.norm_rmvfilt_rmvsamp,
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12,
         filename = "pheatmap_z_bcountdatauq.bouschetmarkers1.norm_rmvfilt_rmvsamp.png")
dev.off()

# Compute distances and hierarchical clustering
bcountdatauq.bouschetmarkers1.norm_rmvfilt <- read.table("bcountdatauq.bouschetmarkers1.norm_rmvfilt.txt", header = T, row.names = 1)
dim(bcountdatauq.bouschetmarkers1.norm_rmvfilt)
head(bcountdatauq.bouschetmarkers1.norm_rmvfilt)
##Basilia Suggestion:Hybrid Type1 category
#Remove samples E14 cell types, ferguson E14 cell types
Tbcountdatauq.bouschetmarkers1.norm_rmvfilt_type1 <- t(bcountdatauq.bouschetmarkers1.norm_rmvfilt[,c(3:12,22:25)])
bcountdatauq.bouschetmarkers1.norm_rmvfilt_type1_scale <- scale(Tbcountdatauq.bouschetmarkers1.norm_rmvfilt_type1)

ddbcountdatauq.bouschetmarkers1.norm_rmvfilt_type1_scale <- dist(as.matrix(bcountdatauq.bouschetmarkers1.norm_rmvfilt_type1_scale), method = "euclidean")
ddbcountdatauq.bouschetmarkers1.norm_rmvfilt_type1_scalehc <- hclust(ddbcountdatauq.bouschetmarkers1.norm_rmvfilt_type1_scale, method = "ward.D")
plot(ddbcountdatauq.bouschetmarkers1.norm_rmvfilt_type1_scalehc)
#save as ddbcountdatauq.bouschetmarkers1.norm_rmvfilt_type1_scalehc.svg


#z-score convert
z_Tbcountdatauq.bouschetmarkers1.norm_rmvfilt_type1= scale(t(bcountdatauq.bouschetmarkers1.norm_rmvfilt[,c(3:12,22:25)]), center = TRUE, scale = TRUE)
z_bcountdatauq.bouschetmarkers1.norm_rmvfilt_type1 <- t(z_Tbcountdatauq.bouschetmarkers1.norm_rmvfilt_type1)

head(z_bcountdatauq.bouschetmarkers1.norm_rmvfilt_type1)
dim(z_bcountdatauq.bouschetmarkers1.norm_rmvfilt_type1)
#filtered samples
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_bcountdatauq.bouschetmarkers1.norm_rmvfilt_type1,
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12,
         filename = "pheatmap_z_bcountdatauq.bouschetmarkers1.norm_rmvfilt_type1.png")
dev.off()

##Basilia Suggestion: type2 category
#keep samples E14 cell types, ferguson E14 cell types and in vitro E14 bouschet
Tbcountdatauq.bouschetmarkers1.norm_rmvfilt_type2 <- t(bcountdatauq.bouschetmarkers1.norm_rmvfilt[,c(1:15,22:25,28:31)])
bcountdatauq.bouschetmarkers1.norm_rmvfilt_type2_scale <- scale(Tbcountdatauq.bouschetmarkers1.norm_rmvfilt_type2)

ddbcountdatauq.bouschetmarkers1.norm_rmvfilt_type2_scale <- dist(as.matrix(bcountdatauq.bouschetmarkers1.norm_rmvfilt_type2_scale), method = "euclidean")
ddbcountdatauq.bouschetmarkers1.norm_rmvfilt_type2_scalehc <- hclust(ddbcountdatauq.bouschetmarkers1.norm_rmvfilt_type2_scale, method = "ward.D")
plot(ddbcountdatauq.bouschetmarkers1.norm_rmvfilt_type2_scalehc)
#save as ddbcountdatauq.bouschetmarkers1.norm_rmvfilt_type2_scalehc.svg
library(ggdendro)
library(ggplot2)
fviz_dend(ddbcountdatauq.bouschetmarkers1.norm_rmvfilt_type2_scalehc, 
          cex = 1, 
          k = 3, 
          k_colors = c("#D95F02", "#D95F02", "#6B8E23"),
          palette= c("#D95F02", "#D95F02", "#6B8E23"),
          horiz = TRUE,
          color_labels_by_k = TRUE,
          type = "rectangle",
          labels_track_height= 25)
#save as fvizddbcountdatauq.bouschetmarkers1.norm_rmvfilt_type2_scalehc.png

#z-score convert
z_Tbcountdatauq.bouschetmarkers1.norm_rmvfilt_type2= scale(t(bcountdatauq.bouschetmarkers1.norm_rmvfilt[,c(1:15,22:25,28:31)]), center = TRUE, scale = TRUE)
z_bcountdatauq.bouschetmarkers1.norm_rmvfilt_type2 <- t(z_Tbcountdatauq.bouschetmarkers1.norm_rmvfilt_type2)

head(z_bcountdatauq.bouschetmarkers1.norm_rmvfilt_type2)
dim(z_bcountdatauq.bouschetmarkers1.norm_rmvfilt_type2)
#filtered samples
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_bcountdatauq.bouschetmarkers1.norm_rmvfilt_type2,
         color = colorRampPalette(c("#2166AC","white","#B2182B"))(length(breaksList1)),
         breaks = breaksList1,
         treeheight_row = 30,
         treeheight_col = 30,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 10, 
         cellheight = 10, filename = "pheatmap_z_bcountdatauq.bouschetmarkers1.norm_rmvfilt_type2.png")
dev.off()



meancheckuq <- bcountdatauq.cellmarkers1.norm_rmvfilt[,c(1:15,22:25,28:31)]
head(meancheckuq)
#transpose
Tmeancheckuq <- t(meancheckuq)
head(Tmeancheckuq)
meanTmeancheckuq <- colMeans(Tmeancheckuq)
sdTmeancheckuq <- colSds(Tmeancheckuq)

TmeancheckuqX <- cbind.data.frame((data.frame(Tmeancheckuq[,1] - meanTmeancheckuq[1])/sdTmeancheckuq[1]),
                                  (data.frame(Tmeancheckuq[,2] - meanTmeancheckuq[2])/sdTmeancheckuq[2]),
                                  (data.frame(Tmeancheckuq[,3] - meanTmeancheckuq[3])/sdTmeancheckuq[3]),
                                  (data.frame(Tmeancheckuq[,4] - meanTmeancheckuq[4])/sdTmeancheckuq[4]),
                                  (data.frame(Tmeancheckuq[,5] - meanTmeancheckuq[5])/sdTmeancheckuq[5]),
                                  (data.frame(Tmeancheckuq[,6] - meanTmeancheckuq[6])/sdTmeancheckuq[6]),
                                  (data.frame(Tmeancheckuq[,7] - meanTmeancheckuq[7])/sdTmeancheckuq[7]),
                                  (data.frame(Tmeancheckuq[,8] - meanTmeancheckuq[8])/sdTmeancheckuq[8]),
                                  (data.frame(Tmeancheckuq[,9] - meanTmeancheckuq[9])/sdTmeancheckuq[9]),
                                  (data.frame(Tmeancheckuq[,10] - meanTmeancheckuq[10])/sdTmeancheckuq[10]),
                                  (data.frame(Tmeancheckuq[,11] - meanTmeancheckuq[11])/sdTmeancheckuq[11]),
                                  (data.frame(Tmeancheckuq[,12] - meanTmeancheckuq[12])/sdTmeancheckuq[12]),
                                  (data.frame(Tmeancheckuq[,13] - meanTmeancheckuq[13])/sdTmeancheckuq[13]),
                                  (data.frame(Tmeancheckuq[,14] - meanTmeancheckuq[14])/sdTmeancheckuq[14]),
                                  (data.frame(Tmeancheckuq[,15] - meanTmeancheckuq[15])/sdTmeancheckuq[15]),
                                  (data.frame(Tmeancheckuq[,16] - meanTmeancheckuq[16])/sdTmeancheckuq[16]),
                                  (data.frame(Tmeancheckuq[,17] - meanTmeancheckuq[17])/sdTmeancheckuq[17]),
                                  (data.frame(Tmeancheckuq[,18] - meanTmeancheckuq[18])/sdTmeancheckuq[18]),
                                  (data.frame(Tmeancheckuq[,19] - meanTmeancheckuq[19])/sdTmeancheckuq[19]),
                                  (data.frame(Tmeancheckuq[,20] - meanTmeancheckuq[20])/sdTmeancheckuq[20]),
                                  (data.frame(Tmeancheckuq[,21] - meanTmeancheckuq[21])/sdTmeancheckuq[21]),
                                  (data.frame(Tmeancheckuq[,22] - meanTmeancheckuq[22])/sdTmeancheckuq[22]),
                                  (data.frame(Tmeancheckuq[,23] - meanTmeancheckuq[23])/sdTmeancheckuq[23]),
                                  (data.frame(Tmeancheckuq[,24] - meanTmeancheckuq[24])/sdTmeancheckuq[24]),
                                  (data.frame(Tmeancheckuq[,25] - meanTmeancheckuq[25])/sdTmeancheckuq[25]),
                                  (data.frame(Tmeancheckuq[,26] - meanTmeancheckuq[26])/sdTmeancheckuq[26]),
                                  (data.frame(Tmeancheckuq[,27] - meanTmeancheckuq[27])/sdTmeancheckuq[27]),
                                  (data.frame(Tmeancheckuq[,28] - meanTmeancheckuq[28])/sdTmeancheckuq[28]),
                                  (data.frame(Tmeancheckuq[,29] - meanTmeancheckuq[29])/sdTmeancheckuq[29]),
                                  (data.frame(Tmeancheckuq[,30] - meanTmeancheckuq[30])/sdTmeancheckuq[30]),
                                  (data.frame(Tmeancheckuq[,31] - meanTmeancheckuq[31])/sdTmeancheckuq[31]),
                                  (data.frame(Tmeancheckuq[,32] - meanTmeancheckuq[32])/sdTmeancheckuq[32]),
                                  (data.frame(Tmeancheckuq[,33] - meanTmeancheckuq[33])/sdTmeancheckuq[33]),
                                  (data.frame(Tmeancheckuq[,34] - meanTmeancheckuq[34])/sdTmeancheckuq[34]),
                                  (data.frame(Tmeancheckuq[,35] - meanTmeancheckuq[35])/sdTmeancheckuq[35]),
                                  (data.frame(Tmeancheckuq[,36] - meanTmeancheckuq[36])/sdTmeancheckuq[36]),
                                  (data.frame(Tmeancheckuq[,37] - meanTmeancheckuq[37])/sdTmeancheckuq[37]),
                                  (data.frame(Tmeancheckuq[,38] - meanTmeancheckuq[38])/sdTmeancheckuq[38]),
                                  (data.frame(Tmeancheckuq[,39] - meanTmeancheckuq[39])/sdTmeancheckuq[39]),
                                  (data.frame(Tmeancheckuq[,40] - meanTmeancheckuq[40])/sdTmeancheckuq[40]),
                                  (data.frame(Tmeancheckuq[,41] - meanTmeancheckuq[41])/sdTmeancheckuq[41]),
                                  (data.frame(Tmeancheckuq[,42] - meanTmeancheckuq[42])/sdTmeancheckuq[42]),
                                  (data.frame(Tmeancheckuq[,43] - meanTmeancheckuq[43])/sdTmeancheckuq[43]),
                                  (data.frame(Tmeancheckuq[,44] - meanTmeancheckuq[44])/sdTmeancheckuq[44]),
                                  (data.frame(Tmeancheckuq[,45] - meanTmeancheckuq[45])/sdTmeancheckuq[45]),
                                  (data.frame(Tmeancheckuq[,46] - meanTmeancheckuq[46])/sdTmeancheckuq[46]),
                                  (data.frame(Tmeancheckuq[,47] - meanTmeancheckuq[47])/sdTmeancheckuq[47]),
                                  (data.frame(Tmeancheckuq[,48] - meanTmeancheckuq[48])/sdTmeancheckuq[48]),
                                  (data.frame(Tmeancheckuq[,49] - meanTmeancheckuq[49])/sdTmeancheckuq[49]),
                                  (data.frame(Tmeancheckuq[,50] - meanTmeancheckuq[50])/sdTmeancheckuq[50]),
                                  (data.frame(Tmeancheckuq[,51] - meanTmeancheckuq[51])/sdTmeancheckuq[51]))

zmeancheckuqX <- t(TmeancheckuqX)
head(zmeancheckuqX,1)
rownames(zmeancheckuqX) <- rownames(meancheckuq)


mediancheckuq <- bcountdatauq.cellmarkers1.norm_rmvfilt[,c(1:15,22:25,28:31)]
head(mediancheckuq)
#transpose
Tmediancheckuq <- t(mediancheckuq)
head(Tmediancheckuq)
medianTmediancheckuq <- colMedians(Tmediancheckuq)
sdTmediancheckuq <- colSds(Tmediancheckuq)

TmediancheckuqX <- cbind.data.frame((data.frame(Tmediancheckuq[,1] - medianTmediancheckuq[1])/sdTmediancheckuq[1]),
                                    (data.frame(Tmediancheckuq[,2] - medianTmediancheckuq[2])/sdTmediancheckuq[2]),
                                    (data.frame(Tmediancheckuq[,3] - medianTmediancheckuq[3])/sdTmediancheckuq[3]),
                                    (data.frame(Tmediancheckuq[,4] - medianTmediancheckuq[4])/sdTmediancheckuq[4]),
                                    (data.frame(Tmediancheckuq[,5] - medianTmediancheckuq[5])/sdTmediancheckuq[5]),
                                    (data.frame(Tmediancheckuq[,6] - medianTmediancheckuq[6])/sdTmediancheckuq[6]),
                                    (data.frame(Tmediancheckuq[,7] - medianTmediancheckuq[7])/sdTmediancheckuq[7]),
                                    (data.frame(Tmediancheckuq[,8] - medianTmediancheckuq[8])/sdTmediancheckuq[8]),
                                    (data.frame(Tmediancheckuq[,9] - medianTmediancheckuq[9])/sdTmediancheckuq[9]),
                                    (data.frame(Tmediancheckuq[,10] - medianTmediancheckuq[10])/sdTmediancheckuq[10]),
                                    (data.frame(Tmediancheckuq[,11] - medianTmediancheckuq[11])/sdTmediancheckuq[11]),
                                    (data.frame(Tmediancheckuq[,12] - medianTmediancheckuq[12])/sdTmediancheckuq[12]),
                                    (data.frame(Tmediancheckuq[,13] - medianTmediancheckuq[13])/sdTmediancheckuq[13]),
                                    (data.frame(Tmediancheckuq[,14] - medianTmediancheckuq[14])/sdTmediancheckuq[14]),
                                    (data.frame(Tmediancheckuq[,15] - medianTmediancheckuq[15])/sdTmediancheckuq[15]),
                                    (data.frame(Tmediancheckuq[,16] - medianTmediancheckuq[16])/sdTmediancheckuq[16]),
                                    (data.frame(Tmediancheckuq[,17] - medianTmediancheckuq[17])/sdTmediancheckuq[17]),
                                    (data.frame(Tmediancheckuq[,18] - medianTmediancheckuq[18])/sdTmediancheckuq[18]),
                                    (data.frame(Tmediancheckuq[,19] - medianTmediancheckuq[19])/sdTmediancheckuq[19]),
                                    (data.frame(Tmediancheckuq[,20] - medianTmediancheckuq[20])/sdTmediancheckuq[20]),
                                    (data.frame(Tmediancheckuq[,21] - medianTmediancheckuq[21])/sdTmediancheckuq[21]),
                                    (data.frame(Tmediancheckuq[,22] - medianTmediancheckuq[22])/sdTmediancheckuq[22]),
                                    (data.frame(Tmediancheckuq[,23] - medianTmediancheckuq[23])/sdTmediancheckuq[23]),
                                    (data.frame(Tmediancheckuq[,24] - medianTmediancheckuq[24])/sdTmediancheckuq[24]),
                                    (data.frame(Tmediancheckuq[,25] - medianTmediancheckuq[25])/sdTmediancheckuq[25]),
                                    (data.frame(Tmediancheckuq[,26] - medianTmediancheckuq[26])/sdTmediancheckuq[26]),
                                    (data.frame(Tmediancheckuq[,27] - medianTmediancheckuq[27])/sdTmediancheckuq[27]),
                                    (data.frame(Tmediancheckuq[,28] - medianTmediancheckuq[28])/sdTmediancheckuq[28]),
                                    (data.frame(Tmediancheckuq[,29] - medianTmediancheckuq[29])/sdTmediancheckuq[29]),
                                    (data.frame(Tmediancheckuq[,30] - medianTmediancheckuq[30])/sdTmediancheckuq[30]),
                                    (data.frame(Tmediancheckuq[,31] - medianTmediancheckuq[31])/sdTmediancheckuq[31]),
                                    (data.frame(Tmediancheckuq[,32] - medianTmediancheckuq[32])/sdTmediancheckuq[32]),
                                    (data.frame(Tmediancheckuq[,33] - medianTmediancheckuq[33])/sdTmediancheckuq[33]),
                                    (data.frame(Tmediancheckuq[,34] - medianTmediancheckuq[34])/sdTmediancheckuq[34]),
                                    (data.frame(Tmediancheckuq[,35] - medianTmediancheckuq[35])/sdTmediancheckuq[35]),
                                    (data.frame(Tmediancheckuq[,36] - medianTmediancheckuq[36])/sdTmediancheckuq[36]),
                                    (data.frame(Tmediancheckuq[,37] - medianTmediancheckuq[37])/sdTmediancheckuq[37]),
                                    (data.frame(Tmediancheckuq[,38] - medianTmediancheckuq[38])/sdTmediancheckuq[38]),
                                    (data.frame(Tmediancheckuq[,39] - medianTmediancheckuq[39])/sdTmediancheckuq[39]),
                                    (data.frame(Tmediancheckuq[,40] - medianTmediancheckuq[40])/sdTmediancheckuq[40]),
                                    (data.frame(Tmediancheckuq[,41] - medianTmediancheckuq[41])/sdTmediancheckuq[41]),
                                    (data.frame(Tmediancheckuq[,42] - medianTmediancheckuq[42])/sdTmediancheckuq[42]),
                                    (data.frame(Tmediancheckuq[,43] - medianTmediancheckuq[43])/sdTmediancheckuq[43]),
                                    (data.frame(Tmediancheckuq[,44] - medianTmediancheckuq[44])/sdTmediancheckuq[44]),
                                    (data.frame(Tmediancheckuq[,45] - medianTmediancheckuq[45])/sdTmediancheckuq[45]),
                                    (data.frame(Tmediancheckuq[,46] - medianTmediancheckuq[46])/sdTmediancheckuq[46]),
                                    (data.frame(Tmediancheckuq[,47] - medianTmediancheckuq[47])/sdTmediancheckuq[47]),
                                    (data.frame(Tmediancheckuq[,48] - medianTmediancheckuq[48])/sdTmediancheckuq[48]),
                                    (data.frame(Tmediancheckuq[,49] - medianTmediancheckuq[49])/sdTmediancheckuq[49]),
                                    (data.frame(Tmediancheckuq[,50] - medianTmediancheckuq[50])/sdTmediancheckuq[50]),
                                    (data.frame(Tmediancheckuq[,51] - medianTmediancheckuq[51])/sdTmediancheckuq[51]))

zmediancheckuqX <- t(TmediancheckuqX)
head(zmediancheckuqX,1)
rownames(zmediancheckuqX) <- rownames(mediancheckuq)
head(zmediancheckuqX,1)
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(zmediancheckuqX,
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 3,
         cellwidth = 20, 
         cellheight = 12,
         filename = "pheatmap_zmediancheckuqX.png")


bcoldata <- read.table("bcoldatasep.txt" , header = TRUE, stringsAsFactors = FALSE)
rownames(bcoldata)  <- c("ourT0WTE14","ourT0ZFP57KOE14","ourT0WTJB1","ourT0ZFP57KOJB1","SRR1409932","SRR1409933","SRR1409934","SRR1409935","SRR1409938","SRR1409939","SRR1409940","SRR1409941","SRR1409944","SRR1409945","SRR1409946","SRR1409950","SRR1409951","SRR1409952","SRR1409953","SRR1409954","SRR1409955","ourT12WT1","ourT12WT2","ourT12ZFP57KO1","ourT12ZFP57KO2","SRR5665889","SRR5665890","SRR8329324","SRR8329325","SRR8329327","SRR8329328")
bcoldata[,1]
rownames(bcoldata)
colnames(bcoldata)
bcoldata = bcoldata[,-1]
bcoldata <- data.frame(bcoldata)
head(bcoldata)
bcoldata <- bcoldata[,c("bcondition","batchnruntype","bcell")]
bcoldata$bcondition <- factor(bcoldata$bcondition)
bcoldata$batchnruntype <- factor(bcoldata$batchnruntype)
bcoldata$bcell <- factor(bcoldata$bcell)

bcoldata1  <- bcoldata[1:6,]
all(rownames(bcoldata1) == colnames(bcountdata1)) #should print TRUE
dds1 <- DESeqDataSetFromMatrix(countData =bcountdata1, colData = bcoldata1, design = ~bcondition)
dds1
keep1 <- rowSums(counts(dds1)) >= 5
dds1 <- dds1[keep1,]
dds1 <- DESeq(dds1)
#Export normalized counts
dds1Normcounts <- counts(dds1, normalized=TRUE)
head(dds1Normcounts)
dim(dds1Normcounts)
write.table(dds1Normcounts, "dds1Normcounts.txt", sep="\t", quote=F, col.names=NA)
plotPCA(as.matrix(dds1Normcounts)) #library(EDASeq)

bcoldata2  <- bcoldata[7:23,]
all(rownames(bcoldata2) == colnames(bcountdata2)) #should print TRUE
dds2 <- DESeqDataSetFromMatrix(countData =bcountdata2, colData = bcoldata2, design = ~bcell)
dds2
keep2 <- rowSums(counts(dds2)) >= 10
dds2 <- dds2[keep2,]
dds2 <- DESeq(dds2)
#Export normalized counts
dds2Normcounts <- counts(dds2, normalized=TRUE)
head(dds2Normcounts)
dim(dds2Normcounts)
write.table(dds2Normcounts, "dds2Normcounts.txt", sep="\t", quote=F, col.names=NA)
plotPCA(as.matrix(dds2Normcounts)) #library(EDASeq)



bcoldata3  <- bcoldata[24:27,]
all(rownames(bcoldata3) == colnames(bcountdata3)) #should print TRUE
dds3 <- DESeqDataSetFromMatrix(countData =bcountdata3, colData = bcoldata3, design = ~bcondition)
dds3
keep3 <- rowSums(counts(dds3)) >= 5
dds3 <- dds3[keep3,]
dds3 <- DESeq(dds3)
#Export normalized counts
dds3Normcounts <- counts(dds3, normalized=TRUE)
head(dds3Normcounts)
dim(dds3Normcounts)
write.table(dds3Normcounts, "dds3Normcounts.txt", sep="\t", quote=F, col.names=NA)
plotPCA(as.matrix(dds3Normcounts)) #library(EDASeq)


bcoldata4  <- bcoldata[28:31,]
all(rownames(bcoldata4) == colnames(bcountdata4)) #should print TRUE
dds4 <- DESeqDataSetFromMatrix(countData =bcountdata4, colData = bcoldata4, design = ~bell)
dds4
keep4 <- rowSums(counts(dds4)) >= 5
dds4 <- dds4[keep4,]
dds4 <- DESeq(dds4)
#Export normalized counts
dds4Normcounts <- counts(dds4, normalized=TRUE)
head(dds4Normcounts)
dim(dds4Normcounts)
write.table(dds4Normcounts, "dds4Normcounts.txt", sep="\t", quote=F, col.names=NA)
plotPCA(as.matrix(dds4Normcounts)) #library(EDASeq)




bcoldata5  <- bcoldata[,32:33]
all(rownames(bcoldata5) == colnames(bcountdata5)) #should print TRUE
dds5 <- DESeqDataSetFromMatrix(countData =bcountdata5, colData = bcoldata5, design = ~bcell)
dds5
keep5 <- rowSums(counts(dds5)) >= 5
dds5 <- dds5[keep5,]
dds5 <- DESeq(dds5)
#Export normalized counts
dds5Normcounts <- counts(dds5, normalized=TRUE)
head(dds5Normcounts)
dim(dds5Normcounts)
write.table(dds5Normcounts, "dds5Normcounts.txt", sep="\t", quote=F, col.names=NA)
plotPCA(as.matrix(dds5Normcounts)) #library(EDASeq)



bcountdata <- cbind.data.frame(bcountdata1, bcountdata2, bcountdata3, bcountdata4, bcountdata5)
library(DESeq2)
bcountdata = as.matrix(bcountdata)
head(bcountdata)
dim(bcountdata)
boxplot(bcountdata, ylim=c(0,100))

bcountdatat0avg <- cbind.data.frame(round(rowMeans(bcountdata[,c(1,3)])),
                                    round(rowMeans(bcountdata[,c(2,4)])),
                                    (bcountdata[,5:31]))

write.table(bcountdatat0avg, "bcountdatat0avg.txt", quote = F, append = F, row.names = T)

head(bcountdatat0avg)                                  
colnames(bcountdatat0avg) <- c("ourT0WT","ourT0ZFP57KO","SRR1409932","SRR1409933","SRR1409934","SRR1409935","SRR1409938","SRR1409939","SRR1409940","SRR1409941","SRR1409944","SRR1409945","SRR1409946","SRR1409950","SRR1409951","SRR1409952","SRR1409953","SRR1409954","SRR1409955","ourT12WT1", "ourT12WT2", "ourT12ZFP57KO1", "ourT12ZFP57KO2","SRR5665889","SRR5665890","SRR8329324","SRR8329325","SRR8329327","SRR8329328")
head(bcountdatat0avg) 
keep <- rowSums(bcountdatat0avg) >= 10
filtbcountdatat0 <- bcountdatat0avg[keep,]
boxplot(log10(filtbcountdatat0), ylim=c(0,10), col = c("blue","blue","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","grey","grey","orange","orange","orange","orange"))
#save as filtbcountdatat0.sep.svg
library(EDASeq)
plotPCA(as.matrix(bcountdata))
plotPCA(as.matrix(bcountdata), labels=F, col =  c("cyan","cyan","red","red","blue","blue","darkgreen","darkgreen","blue","blue","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","orange","orange","orange","orange","orange","orange","navyblue","navyblue","magenta","magenta","cyan","cyan","cyan","cyan","red","red"))
ggsave("PCA_tbcountdata_scaleT.sep.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

#Dendrogram cluster
dim(bcountdata)
head(bcountdata)
# Compute distances and hierarchical clustering
bcountdatatall <- t(bcountdata)
#head(bcountdataavg1)
bcountscale <- scale(bcountdatatall)
#head(bcountavgscale)
ddbcount <- dist(as.matrix(bcountscale), method = "euclidean")
#head(ddbcountavg)
bcounthc <- hclust(ddbcount, method = "ward.D2")
plot(bcounthc)

head(bcountdatat0avg)
#Dendrogram cluster
dim(bcountdatat0avg)
# Compute distances and hierarchical clustering
bcountdatat0avg1 <- t(bcountdatat0avg)
#head(bcountdatat0avg1)
bcountt0scale <- scale(bcountdatat0avg1)
#head(bcountt0scale)
ddbcountt0 <- dist(as.matrix(bcountt0scale), method = "euclidean")
#head(ddbcountt0)
bcountt0hc <- hclust(ddbcountt0, method = "ward.D2")
plot(bcountt0hc)
#save as bcountt0hc.sep.svg

#Create Deseq2 matrix
head(bcountdata)
colnames(bcountdata)
bcoldata <- read.table("bcoldatasep.txt" , header = TRUE, stringsAsFactors = FALSE)
rownames(bcoldata)  <- colnames(bcountdata)
bcoldata[,1]
rownames(bcoldata)
colnames(bcoldata)
bcoldata = bcoldata[,-1]
bcoldata <- data.frame(bcoldata)
head(bcoldata)
bcoldata <- bcoldata[,c("bcondition","batchnruntype","bcell")]
bcoldata$bcondition <- factor(bcoldata$bcondition)
bcoldata$batchnruntype <- factor(bcoldata$batchnruntype)
bcoldata$bcell <- factor(bcoldata$bcell)

all(rownames(bcoldata) == colnames(bcountdata)) #should print TRUE

Bdds <- DESeqDataSetFromMatrix(countData =bcountdata, colData = bcoldata, design = ~batchnruntype + bcondition)
Bdds
#featureData <- data.frame(gene=rownames(bcountdata))
bkeep <- rowSums(counts(Bdds)) >= 100
Bdds <- Bdds[bkeep,]

#View filtered count matrix: View(counts(gdds))
#Apply first method Size factor based of DESeq2
BddsSF <- Bdds
BddsSF <- estimateSizeFactors(BddsSF)
sizeFactors(BddsSF)
Bdds <-DESeq(Bdds)

#Collect most variable genes compare ESC WT n NPC WT
#Extract WTs
head(bcountdata,2)
bcountdataWT <- bcountdata[,c(1:3,7:25,28:31)]
head(bcountdataWT)
bcoldataWT <- bcoldata[c(1:3,7:25,28:31),]
BddsWT <- DESeqDataSetFromMatrix(countData =bcountdataWT, colData = bcoldataWT, design = ~batchnruntype + bcell)
BddsWT
#featureData <- data.frame(gene=rownames(bcountdata))
bkeepWT <- rowSums(counts(BddsWT)) >= 10
BddsWT <- BddsWT[bkeepWT,]
BddsWT <-DESeq(BddsWT)

#DEG
#Normalization is the part of DESeq command: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#See size factors estimated by DESeq function: sizeFactors(dds)
res <- results(BddsWT, contrast=c("bcell","NPC","ESC"))
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
dim(res_sorted)
res_sorted_gene <- data.frame(cbind(chr.pos.1, res_sorted))
head(res_sorted_gene)
dim(res_sorted_gene)
res_sorted_gene_rearranged <- res_sorted_gene[,c(3:5, 1:2, 6:12)]
res_sorted_gene_rearranged <- res_sorted_gene_rearranged[order(res_sorted_gene_rearranged$chr, res_sorted_gene_rearranged$start),]
head(res_sorted_gene_rearranged)
dim(res_sorted_gene_rearranged)
write.table(res_sorted_gene_rearranged, "res_sorted_gene_rearranged.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F, col.names = F)
res_sorted_gene_rearranged_mostvar <- res_sorted_gene_rearranged[order(res_sorted_gene_rearranged$padj),]
head(res_sorted_gene_rearranged_mostvar)
write.table(res_sorted_gene_rearranged_mostvar, "res_sorted_gene_rearranged_mostvar.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F, col.names = F)

Bvsd <- vst(Bdds, blind=FALSE)
head(assay(Bvsd), 3)
sampleDists <- dist(t(assay(Bvsd))) #default is method = "euclidean"
sampleDistshc = hclust(sampleDists, method = "ward.D")
plot(sampleDistshc)
plotPCA(assay(Bvsd), labels=F,col =  c("cyan","cyan","red","red","blue","blue","darkgreen","darkgreen","blue","blue","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","orange","orange","orange","orange","orange","orange","navyblue","navyblue","darkred","darkred","cyan","cyan","cyan","cyan","red","red"))

#With t0 avg n round

bcoldatat0 <- read.table("bcoldatat0.txt" , header = TRUE, stringsAsFactors = FALSE)
rownames(bcoldatat0)  <- colnames(bcountdatat0avg)
bcoldatat0[,1]
rownames(bcoldatat0)
colnames(bcoldatat0)
bcoldatat0 = bcoldatat0[,-1]
bcoldatat0 <- data.frame(bcoldatat0)
head(bcoldatat0)
bcoldatat0 <- bcoldatat0[,c("bcondition","batchnruntype", "bcell")]
bcoldatat0$bcondition <- factor(bcoldatat0$bcondition)
bcoldatat0$batchnruntype <- factor(bcoldatat0$batchnruntype)
bcoldatat0$bcell <- factor(bcoldatat0$bcell)
all(rownames(bcoldatat0) == colnames(bcountdatat0avg)) #should print TRUE
head(bcountdatat0avg)
Bddst0 <- DESeqDataSetFromMatrix(countData =bcountdatat0avg, colData = bcoldatat0, design = ~batchnruntype + bcell)
Bddst0
#featureData <- data.frame(gene=rownames(bcountdata))
bt0keep <- rowSums(counts(Bddst0)) >= 10
Bddst0 <- Bddst0[bt0keep,]

#View filtered count matrix: View(counts(gdds))
#Apply first method Size factor based of DESeq2
Bddst0SF <- Bddst0
Bddst0SF <- estimateSizeFactors(Bddst0SF)
sizeFactors(Bddst0SF)
Bddst0 <- DESeq(Bddst0)
Bvsdt0 <- vst(Bddst0, blind=FALSE)
head(assay(Bvsdt0), 3)
sampleDistst0 <- dist(t(assay(Bvsdt0))) #default is method = "euclidean"
sampleDistst0hc = hclust(sampleDistst0, method = "ward.D")
plot(sampleDistst0hc)
#save as sampleDistst0hc.sep.svg
plotPCA(assay(Bvsdt0), labels=F,col =  c("cyan","red","cyan","red","blue","blue","darkgreen","darkgreen","blue","blue","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","orange","orange","orange","orange","orange","orange","navyblue","navyblue","darkred","darkred","darkgreen","darkgreen","cyan","cyan","red","red"))

#xxxxxxxxx

BddsSFnormcounts <- counts(BddsSF, normalized=TRUE)
head(BddsSFnormcounts)
BddsSFnormcounts <- data.frame(BddsSFnormcounts[order(rownames(BddsSFnormcounts)),])
head(BddsSFnormcounts)
dim(BddsSFnormcounts)

colSums(BddsSFnormcounts)
boxplot(BddsSFnormcounts, ylim=c(0,1))
plotPCA(as.matrix(BddsSFnormcounts))
plotPCA(as.matrix(BddsSFnormcounts), labels=F, col =  c("cyan","cyan","red","red","blue","blue","darkgreen","darkgreen","blue","blue","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","orange","orange","orange","orange","orange","orange","navyblue","navyblue","darkred","darkred","cyan","cyan","cyan","cyan","red","red"))
ggsave("PCA_BddsSFnormcounts_scaleT.sep.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

#Dendrogram cluster
dim(BddsSFnormcounts)
# Compute distances and hierarchical clustering
BddsSFnormcounts1 <- t(BddsSFnormcounts)
BddsSFnormcountsscale <- scale(BddsSFnormcounts1)

ddBddsSFnormcountsscale <- dist(as.matrix(BddsSFnormcountsscale), method = "euclidean")
ddBddsSFnormcountsscalehc <- hclust(ddBddsSFnormcountsscale, method = "ward.D2")
plot(ddBddsSFnormcountsscalehc)


#Average our T0 ESCs and keep others replicates
head(BddsSFnormcounts)
BddsSFt0avgnormcounts <- cbind.data.frame(rowMeans(BddsSFnormcounts[,1:3]),
                                          rowMeans(BddsSFnormcounts[,4:6]),
                                          BddsSFnormcounts[,7:33])
head(BddsSFt0avgnormcounts)
dim(BddsSFt0avgnormcounts)
colnames(BddsSFt0avgnormcounts) <- c("ourT0WT","ourT0ZFP57KO","SRR1409932","SRR1409933","SRR1409934","SRR1409935","SRR1409938","SRR1409939","SRR1409940","SRR1409941","SRR1409944","SRR1409945","SRR1409946","SRR1409950","SRR1409951","SRR1409952","SRR1409953","SRR1409954","SRR1409955","ourT12WT1", "ourT12WT2", "ourT12ZFP57KO1", "ourT12ZFP57KO2","SRR5665889","SRR5665890","SRR8329324","SRR8329325","SRR8329327","SRR8329328")
write.table(BddsSFt0avgnormcounts,"BddsSFt0avgnormcounts.txt", quote = F, append = F, row.names = T)
colSums(BddsSFt0avgnormcounts)
boxplot(BddsSFt0avgnormcounts, ylim=c(0,10))
plotPCA(as.matrix(BddsSFt0avgnormcounts))
plotPCA(as.matrix(BddsSFt0avgnormcounts), labels=F, col =  c("cyan","red","cyan","red","blue","blue","darkgreen","darkgreen","blue","blue","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","orange","orange","orange","orange","orange","orange","navyblue","navyblue","darkred","darkred","cyan","cyan","cyan","cyan","red","red"))
ggsave("PCA_BddsSFt0avgnormcounts_scaleT.sep.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

#Dendrogram cluster
dim(BddsSFt0avgnormcounts)
# Compute distances and hierarchical clustering
BddsSFt0avgnormcounts1 <- t(BddsSFt0avgnormcounts)
BddsSFt0avgnormcountsscale <- scale(BddsSFt0avgnormcounts1)

ddBddsSFt0avgnormcountsscale <- dist(as.matrix(BddsSFt0avgnormcountsscale), method = "euclidean")
ddBddsSFt0avgnormcountsscalehc <- hclust(ddBddsSFt0avgnormcountsscale, method = "ward.D2")
plot(ddBddsSFt0avgnormcountsscalehc)
#save as ddBddsSFt0avgnormcountsscalehc.sep.svg

#with T0 already avg n round 
BddsSFt0normcounts <- counts(Bddst0SF, normalized=TRUE)
head(BddsSFt0normcounts)
BddsSFt0normcounts <- data.frame(BddsSFt0normcounts[order(rownames(BddsSFt0normcounts)),])
head(BddsSFt0normcounts)
dim(BddsSFt0normcounts)
write.table(BddsSFt0normcounts, "BddsSFt0normcounts.txt", quote = F, append = F, row.names = T)
colSums(BddsSFt0normcounts)
boxplot(BddsSFt0normcounts, ylim=c(0,1))

head(BddsSFt0normcounts) 
keept0 <- rowSums(BddsSFt0normcounts) >= 100
filtBddsSFt0normt0 <- BddsSFt0normcounts[keept0,]
boxplot(log10(filtBddsSFt0normt0), ylim=c(0,10), col = c("blue","blue","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","grey","grey","orange","orange","orange","orange"))
#save as filtBddsSFt0normt0.sep.svg

#Dendrogram cluster
dim(BddsSFt0normcounts)
# Compute distances and hierarchical clustering
BddsSFt0normcounts1 <- t(BddsSFt0normcounts)
BddsSFt0normcountsscale <- scale(BddsSFt0normcounts1)

ddBddsSFt0normcountsscale <- dist(as.matrix(BddsSFt0normcountsscale), method = "euclidean")
ddBddsSFt0normcountsscalehc <- hclust(ddBddsSFt0normcountsscale, method = "ward.D2")
plot(ddBddsSFt0normcountsscalehc)
#save as ddBddsSFt0normcountsscalehc.sep.svg

head(BddsSFt0avgnormcounts)
write.table(BddsSFt0avgnormcounts, "BddsSFt0avgnormcounts.txt", quote=F, append=F, row.names = T)
dim(BddsSFt0avgnormcounts)


head(BddsSFt0avgnormcounts)
BddsSFt0avgnormcounts["id"] <- rownames(BddsSFt0avgnormcounts)

#fgrep -f cellmarkers.symbol /home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt -w | grep ps -v  | awk '{print $3"\t"$4"\t"$5"\t"$1"\t"$2}' > cellmarkers.symbol.ens.chr.txt


cellmarkers.cords <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/cellmarkers.symbol.ens.chr_filt.txt", header = F)
head(cellmarkers.cords)
colnames(cellmarkers.cords) <- c("chr", "start", "end", "id", "Gene")
BddsSFt0avgnormcounts.cellmarkers = merge(BddsSFt0avgnormcounts, cellmarkers.cords, by="id", all.x=FALSE)
head(BddsSFt0avgnormcounts.cellmarkers)
dim(BddsSFt0avgnormcounts.cellmarkers)
BddsSFt0avgnormcounts.cellmarkers <- BddsSFt0avgnormcounts.cellmarkers[,c(31:34, 1, 2:30)]
head(BddsSFt0avgnormcounts.cellmarkers,2)
rownames(BddsSFt0avgnormcounts.cellmarkers) <- BddsSFt0avgnormcounts.cellmarkers$Gene
BddsSFt0avgnormcounts.cellmarkers1 <- BddsSFt0avgnormcounts.cellmarkers[,c(6:34)]
head(BddsSFt0avgnormcounts.cellmarkers1,2)
BddsSFt0avgnormcounts.cellmarkers1 <- as.matrix(BddsSFt0avgnormcounts.cellmarkers1)

#Loading control
loadingcontrol.cords <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/loadingcontrol.symbol.ens.chr.txt", header = F)
head(loadingcontrol.cords)
colnames(loadingcontrol.cords) <- c("chr", "start", "end", "id", "Gene")
BddsSFt0avgnormcounts.loadingcontrol = merge(BddsSFt0avgnormcounts, loadingcontrol.cords, by="id", all.x=FALSE)
head(BddsSFt0avgnormcounts.loadingcontrol)
dim(BddsSFt0avgnormcounts.loadingcontrol)
BddsSFt0avgnormcounts.loadingcontrol <- BddsSFt0avgnormcounts.loadingcontrol[,c(31:34, 1, 2:30)]
head(BddsSFt0avgnormcounts.loadingcontrol,2)
rownames(BddsSFt0avgnormcounts.loadingcontrol) <- BddsSFt0avgnormcounts.loadingcontrol$Gene
BddsSFt0avgnormcounts.loadingcontrol1 <- BddsSFt0avgnormcounts.loadingcontrol[,c(6:34)]
head(BddsSFt0avgnormcounts.loadingcontrol1,2)
BddsSFt0avgnormcounts.loadingcontrol1 <- as.matrix(BddsSFt0avgnormcounts.loadingcontrol1)
BddsSFt0avgnormcounts.loadingcontrol1.Actb <- BddsSFt0avgnormcounts.loadingcontrol1[1,]
#Normalize with loading control

BddsSFt0avgnormcounts.cellmarkers1.norm <- cbind.data.frame((BddsSFt0avgnormcounts.cellmarkers1[,1] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[1]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,2] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[2]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,3] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[3]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,4] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[4]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,5] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[5]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,6] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[6]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,7] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[7]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,8] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[8]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,9] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[9]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,10] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[10]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,11] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[11]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,12] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[12]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,13] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[13]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,14] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[14]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,15] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[15]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,16] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[16]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,17] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[17]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,18] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[18]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,19] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[19]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,20] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[20]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,21] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[21]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,22] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[22]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,23] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[23]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,24] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[24]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,25] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[25]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,26] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[26]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,27] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[27]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,28] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[28]),
                                                            (BddsSFt0avgnormcounts.cellmarkers1[,29] / BddsSFt0avgnormcounts.loadingcontrol1.Actb[29]))


head(BddsSFt0avgnormcounts.cellmarkers1.norm)                                                             
colnames(BddsSFt0avgnormcounts.cellmarkers1.norm) <- colnames(BddsSFt0avgnormcounts.cellmarkers1)
head(BddsSFt0avgnormcounts.cellmarkers1.norm, 3)                                                             
z_TBddsSFt0avgnormcounts.cellmarkers1.norm= scale(t(BddsSFt0avgnormcounts.cellmarkers1.norm), center = TRUE, scale = TRUE)
z_BddsSFt0avgnormcounts.cellmarkers1.norm <- t(z_TBddsSFt0avgnormcounts.cellmarkers1.norm)
colfunc <- colorRampPalette(c("navy","white","firebrick3"))

#Dendrogram cluster
dim(BddsSFt0avgnormcounts.cellmarkers1.norm)
# Compute distances and hierarchical clustering
BddsSFt0avgnormcounts.cellmarkers1.norm1 <- t(BddsSFt0avgnormcounts.cellmarkers1.norm)
BddsSFt0avgnormcounts.cellmarkers1.normscale <- scale(BddsSFt0avgnormcounts.cellmarkers1.norm1)

ddBddsSFt0avgnormcounts.cellmarkers1.normscale <- dist(as.matrix(BddsSFt0avgnormcounts.cellmarkers1.normscale), method = "euclidean")
ddBddsSFt0avgnormcounts.cellmarkers1.normscalehc <- hclust(ddBddsSFt0avgnormcounts.cellmarkers1.normscale, method = "ward.D")
plot(ddBddsSFt0avgnormcounts.cellmarkers1.normscalehc)
#save as ddBddsSFt0avgnormcounts.cellmarkers1.normscalehc.sep.svg

svg(filename="heatmap_z_BddsSFt0avgnormcounts.cellmarkers1.norm.sep.svg", width=8, height=18, pointsize=12)
heatmap.2(z_BddsSFt0avgnormcounts.cellmarkers1.norm,trace = "none", col = colfunc,
          density.info=c("none"), dendrogram="both", scale = "none", 
          sepwidth=c(0.001, 0.001), cexRow=2, font=3, cexCol = 0.8, 
          sepcolor="black", margins =c(6,8), srtCol = 45, 
          breaks = seq(-0.1,0.1, length.out = 100), colsep=1:ncol(z_BddsSFt0avgnormcounts.cellmarkers1.norm),
          rowsep=1:nrow(z_BddsSFt0avgnormcounts.cellmarkers1.norm))
dev.off()

library(pheatmap)
library(RColorBrewer)
breaksList = seq(-1, 1)
#pheatmap(data1,treeheight_row = 0, cluster_cols=F, cluster_rows=F, treeheight_col = 0, gaps_col =NULL, gaps_row = NULL, border_color = "black", breaks = breaksList,color= colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)))
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_BddsSFt0avgnormcounts.cellmarkers1.norm,
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 3,
         cellwidth = 20, 
         cellheight = 12, 
         filename = "pheatmap_z_BddsSFt0avgnormcounts.cellmarkers1.norm.sep.png")

#Batch Effect detected#
#Remove batch effect using Combat#
sample_info <- read.table("sampleinfoadd.txt", header = TRUE)
head(sample_info)
sample_info[,1]
rownames(sample_info)=sample_info[,1]
sample_info = sample_info[,-1]
head(sample_info)

phenobatch <- read.table("phenoadd.txt", header = TRUE)
rownames(phenobatch)
phenobatch[,1]
rownames(phenobatch)=phenobatch[,1]
rownames(phenobatch)
colnames(phenobatch)
phenobatch = phenobatch[,-1]
head(phenobatch)
colnames(phenobatch)
dim(phenobatch)
Batch = data.frame(phenobatch$Batch)
head(Batch)
rownames(Batch) <- colnames(BddsSFt0avgnormcounts[,1:29])
head(Batch)
Run = data.frame(phenobatch$Run)
rownames(Run) <- colnames(BddsSFt0avgnormcounts[,1:29])
head(Run)
Modcombat = model.matrix(~1, data=Batch)

#BiocManager::install("sva")
library(sva)
BddsSFt0avgnormcountsfilt <- BddsSFt0avgnormcounts[,1:29]
BddsSFt0avgnormcountsfilt[BddsSFt0avgnormcountsfilt==0] <- NA
head(BddsSFt0avgnormcountsfilt)
BddsSFt0avgnormcountsfilt1 <-BddsSFt0avgnormcountsfilt[complete.cases(BddsSFt0avgnormcountsfilt),]
dim(BddsSFt0avgnormcountsfilt1)
head(BddsSFt0avgnormcountsfilt1)

head(BddsSFt0avgnormcountsfilt1,1)
write.table(BddsSFt0avgnormcountsfilt1, "BddsSFt0avgnormcountsfilt1.txt", quote=F, append = F, row.names = T)
BddsSFt0avgnormcountscor = ComBat(dat=as.matrix(BddsSFt0avgnormcountsfilt1), batch=as.numeric(phenobatch$Batch), mod=Modcombat, par.prior=TRUE, prior.plots=FALSE)
#BddsSFt0avgnormcountscor = ComBat(dat=as.matrix(BddsSFt0avgnormcountscorX), batch=as.numeric(phenobatch$Run), mod=Runcombat, par.prior=TRUE, prior.plots=FALSE)

head(BddsSFt0avgnormcountscor)
dim(BddsSFt0avgnormcountscor)
write.table(BddsSFt0avgnormcountscor, "BddsSFt0avgnormcountscor.txt", quote = F, append= F, row.names = T)
boxplot(BddsSFt0avgnormcountscor, ylim=c(0,2000))
plotPCA(as.matrix(BddsSFt0avgnormcountscor))
plotPCA(as.matrix(BddsSFt0avgnormcountscor), labels=F, col =  c("cyan","red","cyan","red","blue","blue","darkgreen","darkgreen","blue","blue","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","orange","orange","orange","orange","orange","orange","navyblue","navyblue","darkred","darkred","cyan","cyan","cyan","cyan","red","red"))
ggsave("PCA_BddsSFt0avgnormcountscor_scaleT.sep.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

#Dendrogram cluster
BddsSFt0avgnormcountscor <- read.table("BddsSFt0avgnormcountscor.txt", header = T)
dim(BddsSFt0avgnormcountscor)
# Compute distances and hierarchical clustering
BddsSFt0avgnormcountscor1 <- t(BddsSFt0avgnormcountscor)
BddsSFt0avgnormcountscorscale <- scale(BddsSFt0avgnormcountscor1)

ddBddsSFt0avgnormcountscorscale <- dist(as.matrix(BddsSFt0avgnormcountscorscale), method = "euclidean")
ddBddsSFt0avgnormcountscorscalehc <- hclust(ddBddsSFt0avgnormcountscorscale, method = "ward.D2")
plot(ddBddsSFt0avgnormcountscorscalehc)
#save as ddBddsSFt0avgnormcountscorscalehc.sep.svg



#fgrep -f /media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/selected.markers.ens BddsSFt0avgnormcountscor.txt > BddsSFt0avgnormcountscor_selectedmarkers.txt

#Dendrogram cluster
BddsSFt0avgnormcountscor_selectedmarkers <- read.table("BddsSFt0avgnormcountscor_selectedmarkers.txt", header =F)
rownames(BddsSFt0avgnormcountscor_selectedmarkers) <- BddsSFt0avgnormcountscor_selectedmarkers[,1]
BddsSFt0avgnormcountscor_selectedmarkers <- BddsSFt0avgnormcountscor_selectedmarkers[,-1]
dim(BddsSFt0avgnormcountscor_selectedmarkers)
head(BddsSFt0avgnormcountscor_selectedmarkers)
colnames(BddsSFt0avgnormcountscor_selectedmarkers) <- c("ourT0WT","ourT0ZFP57KO","SRR1409932","SRR1409933","SRR1409934","SRR1409935","SRR1409938","SRR1409939","SRR1409940","SRR1409941","SRR1409944","SRR1409945","SRR1409946","SRR1409950","SRR1409951","SRR1409952","SRR1409953","SRR1409954","SRR1409955","ourT12WT1", "ourT12WT2", "ourT12ZFP57KO1", "ourT12ZFP57KO2","SRR5665889","SRR5665890","SRR8329324","SRR8329325","SRR8329327","SRR8329328")
# Compute distances and hierarchical clustering
BddsSFt0avgnormcountscor_selectedmarkers1 <- t(BddsSFt0avgnormcountscor_selectedmarkers)
BddsSFt0avgnormcountscor_selectedmarkersscale <- scale(BddsSFt0avgnormcountscor_selectedmarkers1)

ddBddsSFt0avgnormcountscor_selectedmarkersscale <- dist(as.matrix(BddsSFt0avgnormcountscor_selectedmarkersscale), method = "euclidean")
ddBddsSFt0avgnormcountscor_selectedmarkersscalehc <- hclust(ddBddsSFt0avgnormcountscor_selectedmarkersscale, method = "ward.D2")
plot(ddBddsSFt0avgnormcountscor_selectedmarkersscalehc)
#save as ddBddsSFt0avgnormcountscor_selectedmarkersscalehc.sep.svg

library(sva)
BddsSFt0normcountsfilt <- BddsSFt0normcounts
BddsSFt0normcountsfilt[BddsSFt0normcountsfilt==0] <- NA
head(BddsSFt0normcountsfilt)
BddsSFt0normcountsfilt1 <-BddsSFt0normcountsfilt[complete.cases(BddsSFt0normcountsfilt),]
dim(BddsSFt0normcountsfilt1)
head(BddsSFt0normcountsfilt1)

head(BddsSFt0normcountsfilt1,1)
write.table(BddsSFt0normcountsfilt1, "BddsSFt0normcountsfilt1.txt", quote=F, append = F, row.names = T)
BddsSFt0normcountscor = ComBat(dat=as.matrix(BddsSFt0normcountsfilt1), batch=as.numeric(phenobatch$Batch), mod=Modcombat, par.prior=TRUE, prior.plots=FALSE)
head(BddsSFt0normcountscor)
dim(BddsSFt0normcountscor)
write.table(BddsSFt0normcountscor, "BddsSFt0normcountscor.txt", quote = F, append= F, row.names = T)
boxplot(log10(BddsSFt0normcountscor), ylim=c(0,10), col = c("blue","blue","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","grey","grey","orange","orange","orange","orange"))
#save as BddsSFt0normcountscor.sep.svg
plotPCA(as.matrix(BddsSFt0normcountscor))
plotPCA(as.matrix(BddsSFt0normcountscor), labels=F, col =  c("cyan","red","cyan","red","blue","blue","darkgreen","darkgreen","blue","blue","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","orange","orange","orange","orange","orange","orange","navyblue","navyblue","darkred","darkred","cyan","cyan","cyan","cyan","red","red"))
ggsave("PCA_BddsSFt0normcountscor_scaleT.sep.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

#Dendrogram cluster
BddsSFt0normcountscor <- read.table("BddsSFt0normcountscor.txt", header = T)
dim(BddsSFt0normcountscor)
# Compute distances and hierarchical clustering
BddsSFt0normcountscor1 <- t(BddsSFt0normcountscor)
BddsSFt0normcountscorscale <- scale(BddsSFt0normcountscor1)

ddBddsSFt0normcountscorscale <- dist(as.matrix(BddsSFt0normcountscorscale), method = "euclidean")
ddBddsSFt0normcountscorscalehc <- hclust(ddBddsSFt0normcountscorscale, method = "ward.D2")
plot(ddBddsSFt0normcountscorscalehc)
#save as ddBddsSFt0normcountscorscalehc.sep.svg

Bddst0highfilt <- DESeqDataSetFromMatrix(countData =bcountdatat0avg, colData = bcoldatat0, design = ~batch + bcondition)
Bddst0highfilt
#featureData <- data.frame(gene=rownames(bcountdata))
bt0highfiltkeep <- rowSums(counts(Bddst0highfilt)) >= 10
Bddst0highfilt <- Bddst0highfilt[bt0highfiltkeep,]

#View filtered count matrix: View(counts(gdds))
#Apply first method Size factor based of DESeq2
Bddst0highfiltSF <- Bddst0highfilt
Bddst0highfiltSF <- estimateSizeFactors(Bddst0highfiltSF)
sizeFactors(Bddst0highfiltSF)
Bvsdt0highfilt <- vst(Bddst0highfiltSF, blind=FALSE)
head(assay(Bvsdt0highfilt), 3)
boxplot(assay(Bvsdt0highfilt), ylim=c(0,20), col = c("blue","blue","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","grey","grey","orange","orange","orange","orange"))
#save as boxplot_Bvsdt0highfilt.sep.svg
sampleDistst0highfilt <- dist(t(assay(Bvsdt0highfilt))) #default is method = "euclidean"
sampleDistst0highfilthc = hclust(sampleDistst0highfilt, method = "ward.D")
plot(sampleDistst0highfilthc)
#save as sampleDistst0highfilthc.sep.svg

#No  norm
#View filtered count matrix: View(counts(gdds))
#Apply first method Size factor based of DESeq2
Bddst0No <- Bddst0
sizeFactors(Bddst0No)
Bvsdt0No <- vst(Bddst0No, blind=FALSE)
head(assay(Bvsdt0No), 3)
head(assay(Bvsdt0No), 3)
boxplot(assay(Bvsdt0No), ylim=c(0,20), col = c("blue","blue","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","grey","grey","orange","orange","orange","orange"))
#save as boxplot_Bvsdt0No.sep.svg

sampleDistst0No <- dist(t(assay(Bvsdt0No))) #default is method = "euclidean"
sampleDistst0Nohc = hclust(sampleDistst0No, method = "ward.D")
plot(sampleDistst0Nohc)
#save as sampleDistst0Nohc.sep.svg

#CPM Normalization
head(bcountdata)
dim(bcountdata)
bcountdata <- data.frame(bcountdata)
colSums(bcountdata)
bcountdataCPM  <- cbind.data.frame(rownames(bcountdata),
                                   (bcountdata$ourT0WT1 * 1000000) / 7900780,
                                   (bcountdata$ourT0WT2 * 1000000) / 8592481,
                                   (bcountdata$ourT0WT3 * 1000000) / 24700476,
                                   (bcountdata$ourT0ZFP57KO1 * 1000000) / 10944905,
                                   (bcountdata$ourT0ZFP57KO2 * 1000000) / 11928711,
                                   (bcountdata$ourT0ZFP57KO3 * 1000000) / 33651907,
                                   (bcountdata$SRR1409932 * 1000000) / 24344998,
                                   (bcountdata$SRR1409933 * 1000000) / 30037077,
                                   (bcountdata$SRR1409934 * 1000000) / 23901611,
                                   (bcountdata$SRR1409935 * 1000000) / 25204938,
                                   (bcountdata$SRR1409938 * 1000000) / 24688627,
                                   (bcountdata$SRR1409939 * 1000000) / 24435652,
                                   (bcountdata$SRR1409940 * 1000000) / 29696587,
                                   (bcountdata$SRR1409941 * 1000000) / 38930972,
                                   (bcountdata$SRR1409944 * 1000000) / 34126110,
                                   (bcountdata$SRR1409945 * 1000000) / 33764513,
                                   (bcountdata$SRR1409946 * 1000000) / 28003047,
                                   (bcountdata$SRR1409950 * 1000000) / 30108376,
                                   (bcountdata$SRR1409951 * 1000000) / 33716810,
                                   (bcountdata$SRR1409952 * 1000000) / 31440182,
                                   (bcountdata$SRR1409953 * 1000000) / 26186023,
                                   (bcountdata$SRR1409954 * 1000000) / 28177620,
                                   (bcountdata$SRR1409955 * 1000000) / 35447266,
                                   (bcountdata$ourT12WT1 * 1000000) / 7663520,
                                   (bcountdata$ourT12WT2 * 1000000) / 10225378,
                                   (bcountdata$ourT12ZFP57KO1 * 1000000) / 9545119,
                                   (bcountdata$ourT12ZFP57KO2 * 1000000) / 10337082,
                                   (bcountdata$SRR5665889 * 1000000) / 93566041,
                                   (bcountdata$SRR5665890 * 1000000) / 52170182,
                                   (bcountdata$SRR8329324 * 1000000) / 23915981,
                                   (bcountdata$SRR8329325 * 1000000) / 27452485,
                                   (bcountdata$SRR8329327 * 1000000) / 24423677,
                                   (bcountdata$SRR8329328* 1000000) / 25132776)

#Average our T0 ESCs and keep others replicates
head(bcountdataCPM)
colnames(bcountdataCPM) <- c("ENS","ourT0WTE14","ourT0ZFP57KOE14","ourT0WTJB1","ourT0ZFP57KOJB1","SRR1409932","SRR1409933","SRR1409934","SRR1409935","SRR1409938","SRR1409939","SRR1409940","SRR1409941","SRR1409944","SRR1409945","SRR1409946","SRR1409950","SRR1409951","SRR1409952","SRR1409953","SRR1409954","SRR1409955","ourT12WT1", "ourT12WT2", "ourT12ZFP57KO1", "ourT12ZFP57KO2","SRR5665889","SRR5665890","SRR8329324","SRR8329325","SRR8329327","SRR8329328")
rownames(bcountdataCPM) <- bcountdataCPM[,1]
bcountdataCPM <- bcountdataCPM[,-1]
bcountdataCPM <- data.frame(bcountdataCPM)
head(bcountdataCPM)
bCPMt0avg <- cbind.data.frame(rowMeans(bcountdataCPM[,1:3]),
                              rowMeans(bcountdataCPM[,4:6]),
                              bcountdataCPM[,7:33])
head(bCPMt0avg)
dim(bCPMt0avg)
colnames(bCPMt0avg) <- c("ourT0WT","ourT0ZFP57KO","SRR1409932","SRR1409933","SRR1409934","SRR1409935","SRR1409938","SRR1409939","SRR1409940","SRR1409941","SRR1409944","SRR1409945","SRR1409946","SRR1409950","SRR1409951","SRR1409952","SRR1409953","SRR1409954","SRR1409955","ourT12WT1", "ourT12WT2", "ourT12ZFP57KO1", "ourT12ZFP57KO2","SRR5665889","SRR5665890","SRR8329324","SRR8329325","SRR8329327","SRR8329328")
colSums(bCPMt0avg)
head(bCPMt0avg)
boxplot(bCPMt0avg, ylim=c(0,10000))
plotPCA(as.matrix(bCPMt0avg))
plotPCA(as.matrix(bCPMt0avg), labels=F, col =  c("cyan","red","cyan","red","blue","blue","darkgreen","darkgreen","blue","blue","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","orange","orange","orange","orange","orange","orange","navyblue","navyblue","darkred","darkred","cyan","cyan","cyan","cyan","red","red"))
ggsave("PCA_bCPMt0avg_scaleT.sep.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

#Dendrogram cluster
dim(bCPMt0avg)
# Compute distances and hierarchical clustering
bCPMt0avg1 <- t(bCPMt0avg)
bCPMt0avgscale <- scale(bCPMt0avg1)

ddbCPMt0avgscale <- dist(as.matrix(bCPMt0avgscale), method = "euclidean")
ddbCPMt0avgscalehc <- hclust(ddbCPMt0avgscale, method = "ward.D")
plot(ddbCPMt0avgscalehc)
#save as ddbCPMt0avgscalehc.sep.svg

head(bCPMt0avg)
write.table(bCPMt0avg, "bCPMt0avg.txt", quote=F, append=F, row.names = T)

#Batch Effect detected#
#Remove batch effect using Combat#
samplesinfoadd <- read.table("sampleinfoadd.txt", header = TRUE)
head(samplesinfoadd)
samplesinfoadd[,1]
rownames(samplesinfoadd)=samplesinfoadd[,1]
samplesinfoadd = samplesinfoadd[,-1]
head(samplesinfoadd)

phenoadd <- read.table("phenoadd.txt", header = TRUE)
rownames(phenoadd)
phenoadd[,1]
rownames(phenoadd)=phenoadd[,1]
rownames(phenoadd)
colnames(phenoadd)
phenoadd = phenoadd[,-1]
head(phenoadd)
colnames(phenoadd)
dim(phenoadd)
Batch = data.frame(phenoadd$Batch)
head(Batch)
rownames(Batch) <- colnames(bCPMt0avg)
head(Batch)
Modcombatadd = model.matrix(~1, data=Batch)

#BiocManager::install("sva")
library(sva)
bCPMt0avgfilt <- bCPMt0avg
bCPMt0avgfilt[bCPMt0avgfilt==0] <- NA
head(bCPMt0avgfilt)
bCPMt0avgfilt1 <-bCPMt0avgfilt[complete.cases(bCPMt0avgfilt),]
dim(bCPMt0avgfilt1)
head(bCPMt0avgfilt1)

head(bCPMt0avgfilt1,1)
write.table(bCPMt0avgfilt1, "bCPMt0avgfilt1.txt", quote=F, append = F, row.names = T)
bCPMt0avgcor = ComBat(dat=as.matrix(bCPMt0avgfilt1), batch=as.numeric(phenoadd$Batch), mod=Modcombatadd, par.prior=TRUE, prior.plots=FALSE)
head(bCPMt0avgcor)
dim(bCPMt0avgcor)
write.table(bCPMt0avgcor, "bCPMt0avgcor.txt", quote = F, append= F, row.names = T)
boxplot(bCPMt0avgcor, ylim=c(0,200))
plotPCA(as.matrix(bCPMt0avgcor))
plotPCA(as.matrix(bCPMt0avgcor), labels=F, col =  c("cyan","red","cyan","red","blue","blue","darkgreen","darkgreen","blue","blue","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","orange","orange","orange","orange","orange","orange","navyblue","navyblue","darkred","darkred","cyan","cyan","cyan","cyan","red","red"))
ggsave("PCA_bCPMt0avgcor_scaleT.sep.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

#Dendrogram cluster
bCPMt0avgcor <- read.table("bCPMt0avgcor.txt", header = T)
dim(bCPMt0avgcor)
# Compute distances and hierarchical clustering
bCPMt0avgcor1 <- t(bCPMt0avgcor)
bCPMt0avgcorscale <- scale(bCPMt0avgcor1)

ddbCPMt0avgcorscale <- dist(as.matrix(bCPMt0avgcorscale), method = "euclidean")
ddbCPMt0avgcorscalehc <- hclust(ddbCPMt0avgcorscale, method = "ward.D2")
plot(ddbCPMt0avgcorscalehc)
#save as ddbCPMt0avgcorscalehc.sep.svg

head(bCPMt0avg,3)

bCPMt0avg["id"] <- rownames(bCPMt0avg)
cellmarkers.cords <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/cellmarkers.symbol.ens.chr_filt.txt", header = F)
head(cellmarkers.cords)
colnames(cellmarkers.cords) <- c("chr", "start", "end", "id", "Gene")
bCPMt0avg.cellmarkers = merge(bCPMt0avg, cellmarkers.cords, by="id", all.x=FALSE)
head(bCPMt0avg.cellmarkers)
dim(bCPMt0avg.cellmarkers)
bCPMt0avg.cellmarkers <- bCPMt0avg.cellmarkers[,c(31:34, 1, 2:30)]
head(bCPMt0avg.cellmarkers,2)
rownames(bCPMt0avg.cellmarkers) <- bCPMt0avg.cellmarkers$Gene
bCPMt0avg.cellmarkers1 <- bCPMt0avg.cellmarkers[,c(6:34)]
head(bCPMt0avg.cellmarkers1,2)
bCPMt0avg.cellmarkers1 <- as.matrix(bCPMt0avg.cellmarkers1)

#Loading control
loadingcontrol.cords <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/loadingcontrol.symbol.ens.chr.txt", header = F)
head(loadingcontrol.cords)
colnames(loadingcontrol.cords) <- c("chr", "start", "end", "id", "Gene")
bCPMt0avg.loadingcontrol = merge(bCPMt0avg, loadingcontrol.cords, by="id", all.x=FALSE)
head(bCPMt0avg.loadingcontrol)
dim(bCPMt0avg.loadingcontrol)
bCPMt0avg.loadingcontrol <- bCPMt0avg.loadingcontrol[,c(31:34, 1, 2:30)]
head(bCPMt0avg.loadingcontrol,2)
rownames(bCPMt0avg.loadingcontrol) <- bCPMt0avg.loadingcontrol$Gene
bCPMt0avg.loadingcontrol1 <- bCPMt0avg.loadingcontrol[,c(6:34)]
head(bCPMt0avg.loadingcontrol1,2)
bCPMt0avg.loadingcontrol1 <- as.matrix(bCPMt0avg.loadingcontrol1)
bCPMt0avg.loadingcontrol1.Actb <- bCPMt0avg.loadingcontrol1[1,]
#Normalize with loading control

bCPMt0avg.cellmarkers1.norm <- cbind.data.frame((bCPMt0avg.cellmarkers1[,1] / bCPMt0avg.loadingcontrol1.Actb[1]),
                                                (bCPMt0avg.cellmarkers1[,2] / bCPMt0avg.loadingcontrol1.Actb[2]),
                                                (bCPMt0avg.cellmarkers1[,3] / bCPMt0avg.loadingcontrol1.Actb[3]),
                                                (bCPMt0avg.cellmarkers1[,4] / bCPMt0avg.loadingcontrol1.Actb[4]),
                                                (bCPMt0avg.cellmarkers1[,5] / bCPMt0avg.loadingcontrol1.Actb[5]),
                                                (bCPMt0avg.cellmarkers1[,6] / bCPMt0avg.loadingcontrol1.Actb[6]),
                                                (bCPMt0avg.cellmarkers1[,7] / bCPMt0avg.loadingcontrol1.Actb[7]),
                                                (bCPMt0avg.cellmarkers1[,8] / bCPMt0avg.loadingcontrol1.Actb[8]),
                                                (bCPMt0avg.cellmarkers1[,9] / bCPMt0avg.loadingcontrol1.Actb[9]),
                                                (bCPMt0avg.cellmarkers1[,10] / bCPMt0avg.loadingcontrol1.Actb[10]),
                                                (bCPMt0avg.cellmarkers1[,11] / bCPMt0avg.loadingcontrol1.Actb[11]),
                                                (bCPMt0avg.cellmarkers1[,12] / bCPMt0avg.loadingcontrol1.Actb[12]),
                                                (bCPMt0avg.cellmarkers1[,13] / bCPMt0avg.loadingcontrol1.Actb[13]),
                                                (bCPMt0avg.cellmarkers1[,14] / bCPMt0avg.loadingcontrol1.Actb[14]),
                                                (bCPMt0avg.cellmarkers1[,15] / bCPMt0avg.loadingcontrol1.Actb[15]),
                                                (bCPMt0avg.cellmarkers1[,16] / bCPMt0avg.loadingcontrol1.Actb[16]),
                                                (bCPMt0avg.cellmarkers1[,17] / bCPMt0avg.loadingcontrol1.Actb[17]),
                                                (bCPMt0avg.cellmarkers1[,18] / bCPMt0avg.loadingcontrol1.Actb[18]),
                                                (bCPMt0avg.cellmarkers1[,19] / bCPMt0avg.loadingcontrol1.Actb[19]),
                                                (bCPMt0avg.cellmarkers1[,20] / bCPMt0avg.loadingcontrol1.Actb[20]),
                                                (bCPMt0avg.cellmarkers1[,21] / bCPMt0avg.loadingcontrol1.Actb[21]),
                                                (bCPMt0avg.cellmarkers1[,22] / bCPMt0avg.loadingcontrol1.Actb[22]),
                                                (bCPMt0avg.cellmarkers1[,23] / bCPMt0avg.loadingcontrol1.Actb[23]),
                                                (bCPMt0avg.cellmarkers1[,24] / bCPMt0avg.loadingcontrol1.Actb[24]),
                                                (bCPMt0avg.cellmarkers1[,25] / bCPMt0avg.loadingcontrol1.Actb[25]),
                                                (bCPMt0avg.cellmarkers1[,26] / bCPMt0avg.loadingcontrol1.Actb[26]),
                                                (bCPMt0avg.cellmarkers1[,27] / bCPMt0avg.loadingcontrol1.Actb[27]),
                                                (bCPMt0avg.cellmarkers1[,28] / bCPMt0avg.loadingcontrol1.Actb[28]),
                                                (bCPMt0avg.cellmarkers1[,29] / bCPMt0avg.loadingcontrol1.Actb[29]))

colnames(bCPMt0avg.cellmarkers1.norm) <- c("ourT0WT","ourT0ZFP57KO","SRR1409932","SRR1409933","SRR1409934","SRR1409935","SRR1409938","SRR1409939","SRR1409940","SRR1409941","SRR1409944","SRR1409945","SRR1409946","SRR1409950","SRR1409951","SRR1409952","SRR1409953","SRR1409954","SRR1409955","ourT12WT1", "ourT12WT2", "ourT12ZFP57KO1", "ourT12ZFP57KO2","SRR5665889","SRR5665890","SRR8329324","SRR8329325","SRR8329327","SRR8329328")
write.table(bCPMt0avg.cellmarkers1.norm, "bCPMt0avg.cellmarkers1.norm.txt", row.names = T, quote=F, append=F)
boxplot(log10(bCPMt0avg.cellmarkers1.norm), ylim=c(0,10), col = c("blue","blue","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","grey","grey","orange","orange","orange","orange"))
#save as bCPMt0avg.cellmarkers1.norm.sep.svg
library(EDASeq)
plotPCA(as.matrix(bCPMt0avg.cellmarkers1.norm))
plotPCA(as.matrix(bCPMt0avg.cellmarkers1.norm), labels=F, col =  c("cyan","cyan","cyan","red","red","red","blue","blue","darkgreen","darkgreen","blue","blue","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","orange","orange","orange","orange","orange","orange","navyblue","navyblue","darkred","darkred","cyan","cyan","cyan","cyan","red","red"))
ggsave("PCA_tbCPMt0avg.cellmarkers1.norm_scaleT.sep.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)


head(bCPMt0avg.cellmarkers1.norm)                                                             
colnames(bCPMt0avg.cellmarkers1.norm) <- colnames(bCPMt0avg.cellmarkers1)
head(bCPMt0avg.cellmarkers1.norm, 2)                                                             
z_TbCPMt0avg.cellmarkers1.norm= scale(t(bCPMt0avg.cellmarkers1.norm), center = TRUE, scale = TRUE)
z_bCPMt0avg.cellmarkers1.norm <- t(z_TbCPMt0avg.cellmarkers1.norm)
colfunc <- colorRampPalette(c("navy","white","firebrick3"))


#Dendrogram cluster
dim(bCPMt0avg.cellmarkers1.norm)
# Compute distances and hierarchical clustering
bCPMt0avg.cellmarkers1.norm1 <- t(bCPMt0avg.cellmarkers1.norm[,c(1:23,26:29)])
bCPMt0avg.cellmarkers1.normscale <- scale(bCPMt0avg.cellmarkers1.norm1)

ddbCPMt0avg.cellmarkers1.normscale <- dist(as.matrix(bCPMt0avg.cellmarkers1.normscale), method = "euclidean")
ddbCPMt0avg.cellmarkers1.normscalehc <- hclust(ddbCPMt0avg.cellmarkers1.normscale, method = "ward.D")
plot(ddbCPMt0avg.cellmarkers1.normscalehc)
#save as ddbCPMt0avg.cellmarkers1.normscalehc.sep.svg

svg(filename="heatmap_z_bCPMt0avg.cellmarkers1.norm.sep.svg", width=8, height=18, pointsize=12)
heatmap.2(z_bCPMt0avg.cellmarkers1.norm,trace = "none", col = colfunc,
          density.info=c("none"), dendrogram="both", scale = "none", 
          sepwidth=c(0.001, 0.001), cexRow=2, font=3, cexCol = 0.8, 
          sepcolor="black", margins =c(6,8), srtCol = 45, 
          breaks = seq(-0.1,0.1, length.out = 100), colsep=1:ncol(z_bCPMt0avg.cellmarkers1.norm),
          rowsep=1:nrow(z_bCPMt0avg.cellmarkers1.norm))
dev.off()

library(pheatmap)
library(RColorBrewer)
breaksList = seq(-1, 1)
#pheatmap(data1,treeheight_row = 0, cluster_cols=F, cluster_rows=F, treeheight_col = 0, gaps_col =NULL, gaps_row = NULL, border_color = "black", breaks = breaksList,color= colorRampPalette(c("navy", "white", "firebrick3"))(length(breaksList)))
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_bCPMt0avg.cellmarkers1.norm[,c(1:23,26:29)],
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 1,
         cellwidth = 20, 
         cellheight = 12, 
         filename = "pheatmap_z_bCPMt0avg.cellmarkers1.norm.sep.png")

pheatmap(bCPMt0avg.cellmarkers1.norm[,c(1:23,26:29)],
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12)

#FPKM edgeR
(base) ankitv@ankitv-lab-riccio:/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12$ awk '{print $1}' cellmarkerslist.sorted.txt > cellmarkers.symbol 
(base) ankitv@ankitv-lab-riccio:/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12$ fgrep -f cellmarkers.symbol /home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt -w | grep ps -v  | awk '{print $3"\t"$4"\t"$5"\t"$1"\t"$2}' > cellmarkers.symbol.ens.chr.txt

#This gene length is already given by featurecounts , so use it.
head(bcountdata)
dim(bcountdata)
countforfpkm <- bcountdata
head(countforfpkm)
countforfpkm <- countforfpkm[order(rownames(countforfpkm)),]
head(countforfpkm)
countforfpkm_filt <- countforfpkm[which(rowSums(countforfpkm[,1:31]) >=0),] #Dont filter for now use all the genes
dim(countforfpkm_filt)
head(countforfpkm_filt)
countforfpkm_filt1 <- countforfpkm_filt[,1:31]
head(countforfpkm_filt1)

length <- read.table("Bulk_ourT0ej_star-featureCounts_GRCm38.mm10.txt", header=TRUE, row.names=1)
# Remove first four columns (chr, start, end, strand)
length <- length[ ,5:ncol(length)]
head(length)
gene.length <- length$Length


#Use edgeR for FPKM calculation, the method example is (readcount * 10^9)/(Depth * readlength * normfactor) eg. (1701 * 10^9) / (7661385 * 3262 * 1.0577571) = 64.34682
#edgeR rpkm()
library(edgeR)
dgeforfpkm <- DGEList(counts=countforfpkm_filt1,genes=data.frame(Length=gene.length), group = NULL)
dgeforfpkm <- calcNormFactors(dgeforfpkm)
dgeFPKM <- rpkm(dgeforfpkm, dgeforfpkm$genes$Length)
head(dgeFPKM,2)
write.table(dgeFPKM, "dgeFPKMsep.txt", sep="\t", quote = FALSE, append = FALSE)
dgeFPKM <- data.frame(dgeFPKM)
head(dgeFPKM,3)
dim(dgeFPKM)
dgeFPKM["id"] <- rownames(dgeFPKM)

#sort -k1,1 -u cellmarkerslist.txt | sort -k2,2 > cellmarkerslist.sorted.txt
cellmarkerslist.sorted <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/cellmarkerslist.sorted.txt", header = F)
head(cellmarkerslist.sorted)
colnames(cellmarkerslist.sorted) <- c("Gene", "markers")
dim(cellmarkerslist.sorted) #print 94 genes
ens_gene_names_chrpos_dedup_M20 <- read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt")
head(ens_gene_names_chrpos_dedup_M20)
ens_gene_names_chrpos_dedup_M20 <- ens_gene_names_chrpos_dedup_M20[,c(3,4,5,1,2)]
head(ens_gene_names_chrpos_dedup_M20)
colnames(ens_gene_names_chrpos_dedup_M20) <- c("chr", "start", "end", "id", "Gene")

cellmarkers = merge(cellmarkerslist.sorted, ens_gene_names_chrpos_dedup_M20, by="Gene", all.x=FALSE)
head(cellmarkers)
dim(cellmarkers) #print all 94 genes assigned
cellmarkersre <- cellmarkers[order(cellmarkers$markers),]
head(cellmarkersre)
dim(cellmarkersre)
#cellmarkers.cords <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/cellmarkers.symbol.ens.chr.txt", header = F)
#head(cellmarkers.cords)
#dim(cellmarkers.cords)
#colnames(cellmarkers.cords) <- c("chr", "start", "end", "id", "Gene")

dgeFPKM.cellmarker = merge(dgeFPKM, cellmarkersre, by="id", all.x=FALSE)
head(dgeFPKM.cellmarker)
dim(dgeFPKM.cellmarker) #Print all 94 genes overlapped and assigned
dgeFPKM.cellmarkers <- dgeFPKM.cellmarker[order(dgeFPKM.cellmarker$markers),]
dim(dgeFPKM.cellmarkers)
head(dgeFPKM.cellmarkers)
dgeFPKM.cellmarkersre <- dgeFPKM.cellmarkers[,c(33, 2:32)]
head(dgeFPKM.cellmarkersre,2)
rownames(dgeFPKM.cellmarkersre) <- dgeFPKM.cellmarkersre[,1]
dgeFPKM.cellmarkers1 <- dgeFPKM.cellmarkersre[,c(2:32)]
head(dgeFPKM.cellmarkers1,2)
dgeFPKM.cellmarkers1 <- as.matrix(dgeFPKM.cellmarkers1)

#Loading control
loadingcontrol.cords <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/loadingcontrol.symbol.ens.chr.txt", header = F)
head(loadingcontrol.cords)
colnames(loadingcontrol.cords) <- c("chr", "start", "end", "id", "Gene")
dgeFPKM.loadingcontrol = merge(dgeFPKM, loadingcontrol.cords, by="id", all.x=FALSE)
head(dgeFPKM.loadingcontrol)
dim(dgeFPKM.loadingcontrol)
dgeFPKM.loadingcontrol <- dgeFPKM.loadingcontrol[,c(33, 2:32)]
head(dgeFPKM.loadingcontrol,2)
rownames(dgeFPKM.loadingcontrol) <- dgeFPKM.loadingcontrol$Gene
dgeFPKM.loadingcontrol1 <- dgeFPKM.loadingcontrol[,c(2:32)]
head(dgeFPKM.loadingcontrol1,2)
dgeFPKM.loadingcontrol1 <- as.matrix(dgeFPKM.loadingcontrol1)
dgeFPKM.loadingcontrol1.Actb <- dgeFPKM.loadingcontrol1[1,]
#Normalize with loading control
dim(dgeFPKM.cellmarkers1)

dgeFPKM.cellmarkers1.norm <- cbind.data.frame((dgeFPKM.cellmarkers1[,1] / dgeFPKM.loadingcontrol1.Actb[1]),
                                              (dgeFPKM.cellmarkers1[,2] / dgeFPKM.loadingcontrol1.Actb[2]),
                                              (dgeFPKM.cellmarkers1[,3] / dgeFPKM.loadingcontrol1.Actb[3]),
                                              (dgeFPKM.cellmarkers1[,4] / dgeFPKM.loadingcontrol1.Actb[4]),
                                              (dgeFPKM.cellmarkers1[,5] / dgeFPKM.loadingcontrol1.Actb[5]),
                                              (dgeFPKM.cellmarkers1[,6] / dgeFPKM.loadingcontrol1.Actb[6]),
                                              (dgeFPKM.cellmarkers1[,7] / dgeFPKM.loadingcontrol1.Actb[7]),
                                              (dgeFPKM.cellmarkers1[,8] / dgeFPKM.loadingcontrol1.Actb[8]),
                                              (dgeFPKM.cellmarkers1[,9] / dgeFPKM.loadingcontrol1.Actb[9]),
                                              (dgeFPKM.cellmarkers1[,10] / dgeFPKM.loadingcontrol1.Actb[10]),
                                              (dgeFPKM.cellmarkers1[,11] / dgeFPKM.loadingcontrol1.Actb[11]),
                                              (dgeFPKM.cellmarkers1[,12] / dgeFPKM.loadingcontrol1.Actb[12]),
                                              (dgeFPKM.cellmarkers1[,13] / dgeFPKM.loadingcontrol1.Actb[13]),
                                              (dgeFPKM.cellmarkers1[,14] / dgeFPKM.loadingcontrol1.Actb[14]),
                                              (dgeFPKM.cellmarkers1[,15] / dgeFPKM.loadingcontrol1.Actb[15]),
                                              (dgeFPKM.cellmarkers1[,16] / dgeFPKM.loadingcontrol1.Actb[16]),
                                              (dgeFPKM.cellmarkers1[,17] / dgeFPKM.loadingcontrol1.Actb[17]),
                                              (dgeFPKM.cellmarkers1[,18] / dgeFPKM.loadingcontrol1.Actb[18]),
                                              (dgeFPKM.cellmarkers1[,19] / dgeFPKM.loadingcontrol1.Actb[19]),
                                              (dgeFPKM.cellmarkers1[,20] / dgeFPKM.loadingcontrol1.Actb[20]),
                                              (dgeFPKM.cellmarkers1[,21] / dgeFPKM.loadingcontrol1.Actb[21]),
                                              (dgeFPKM.cellmarkers1[,22] / dgeFPKM.loadingcontrol1.Actb[22]),
                                              (dgeFPKM.cellmarkers1[,23] / dgeFPKM.loadingcontrol1.Actb[23]),
                                              (dgeFPKM.cellmarkers1[,24] / dgeFPKM.loadingcontrol1.Actb[24]),
                                              (dgeFPKM.cellmarkers1[,25] / dgeFPKM.loadingcontrol1.Actb[25]),
                                              (dgeFPKM.cellmarkers1[,26] / dgeFPKM.loadingcontrol1.Actb[26]),
                                              (dgeFPKM.cellmarkers1[,27] / dgeFPKM.loadingcontrol1.Actb[27]),
                                              (dgeFPKM.cellmarkers1[,28] / dgeFPKM.loadingcontrol1.Actb[28]),
                                              (dgeFPKM.cellmarkers1[,29] / dgeFPKM.loadingcontrol1.Actb[29]),
                                              (dgeFPKM.cellmarkers1[,30] / dgeFPKM.loadingcontrol1.Actb[30]),
                                              (dgeFPKM.cellmarkers1[,31] / dgeFPKM.loadingcontrol1.Actb[31]))
colnames(dgeFPKM.cellmarkers1.norm) <- c("ourT0WTE14","ourT0ZFP57KOE14","ourT0WTJB1","ourT0ZFP57KOJB1","SRR1409932","SRR1409933","SRR1409934","SRR1409935","SRR1409938","SRR1409939","SRR1409940","SRR1409941","SRR1409944","SRR1409945","SRR1409946","SRR1409950","SRR1409951","SRR1409952","SRR1409953","SRR1409954","SRR1409955","ourT12WT1", "ourT12WT2", "ourT12ZFP57KO1", "ourT12ZFP57KO2","SRR5665889","SRR5665890","SRR8329324","SRR8329325","SRR8329327","SRR8329328")

head(dgeFPKM.cellmarkers1.norm)
boxplot(log10(dgeFPKM.cellmarkers1.norm), ylim=c(0,10), col = c("blue","blue","blue","blue","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","grey","grey","orange","orange","orange","orange"))
#save as dgeFPKM.cellmarkers1.norm.sep.svg
library(EDASeq)
plotPCA(as.matrix(dgeFPKM.cellmarkers1.norm))
plotPCA(as.matrix(dgeFPKM.cellmarkers1.norm), labels=F, col =  c("cyan","red","cyan","red","blue","blue","darkgreen","darkgreen","blue","blue","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","seagreen","seagreen","seagreen","seagreen","seagreen","seagreen","green","green","magenta","magenta","darkgreen","darkgreen","navy","navy","red","red"))
ggsave("PCA_tdgeFPKM.cellmarkers1.norm_scaleT.sep.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

write.table(dgeFPKM.cellmarkers1.norm, "dgeFPKM.cellmarkers1.norm.sep.txt", row.names = T, quote = F, append=F)
#Dendrogram cluster
dim(dgeFPKM.cellmarkers1.norm)  #94 31
# Compute distances and hierarchical clustering
dgeFPKM.cellmarkers1.norm1 <- t(dgeFPKM.cellmarkers1.norm)
dgeFPKM.cellmarkers1.normscale <- scale(dgeFPKM.cellmarkers1.norm1)
dddgeFPKM.cellmarkers1.normscale <- dist(as.matrix(dgeFPKM.cellmarkers1.normscale), method = "euclidean")
dddgeFPKM.cellmarkers1.normscalehc <- hclust(dddgeFPKM.cellmarkers1.normscale, method = "ward.D")
plot(dddgeFPKM.cellmarkers1.normscalehc)
#save as dddgeFPKM.cellmarkers1.normscalehc.sep.svg
head(dgeFPKM.cellmarkers1.norm,1)                                                             
colnames(dgeFPKM.cellmarkers1.norm) <- colnames(dgeFPKM.cellmarkers1)
head(dgeFPKM.cellmarkers1.norm, 2)                                                             
z_TdgeFPKM.cellmarkers1.norm= scale(t(dgeFPKM.cellmarkers1.norm), center = TRUE, scale = TRUE)
z_dgeFPKM.cellmarkers1.norm <- t(z_TdgeFPKM.cellmarkers1.norm)
colfunc <- colorRampPalette(c("navy","white","firebrick3"))

dim(z_dgeFPKM.cellmarkers1.norm)  #94 31
library(pheatmap)
library(RColorBrewer)
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_dgeFPKM.cellmarkers1.norm[,c(1:31)],
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12, 
         filename = "pheatmap_z_dgeFPKM.cellmarkers1.norm.sep.png")

dev.off()
#Dendrogram cluster
dim(dgeFPKM.cellmarkers1.norm)
head(dgeFPKM.cellmarkers1.norm)
#Filter samples and filter genes
#fgrep -f ./../../2019/time12/time0-tim12/removegenescellmarkers.txt dgeFPKM.cellmarkers1.norm.sep.txt  -v > dgeFPKM.cellmarkers1.norm_rmvfilt.txt
# Compute distances and hierarchical clustering
dgeFPKM.cellmarkers1.norm_rmvfilt <- read.table("dgeFPKM.cellmarkers1.norm_rmvfilt.txt", header = T, row.names = 1)
dim(dgeFPKM.cellmarkers1.norm_rmvfilt)
head(dgeFPKM.cellmarkers1.norm_rmvfilt)
#Remove samples E14 cell types
TdgeFPKM.cellmarkers1.norm_rmvfilt_rmvsamp <- t(dgeFPKM.cellmarkers1.norm_rmvfilt[,c(3:12,22:25,28:31)])
dgeFPKM.cellmarkers1.norm_rmvfilt_rmvsamp_scale <- scale(TdgeFPKM.cellmarkers1.norm_rmvfilt_rmvsamp)

dddgeFPKM.cellmarkers1.norm_rmvfilt_rmvsamp_scale <- dist(as.matrix(dgeFPKM.cellmarkers1.norm_rmvfilt_rmvsamp_scale), method = "euclidean")
dddgeFPKM.cellmarkers1.norm_rmvfilt_rmvsamp_scalehc <- hclust(dddgeFPKM.cellmarkers1.norm_rmvfilt_rmvsamp_scale, method = "ward.D")
plot(dddgeFPKM.cellmarkers1.norm_rmvfilt_rmvsamp_scalehc)
#save as dddgeFPKM.cellmarkers1.norm_rmvfilt_rmvsamp_scalehc.svg


#z-score convert
z_TdgeFPKM.cellmarkers1.norm_rmvfilt_rmvsamp= scale(t(dgeFPKM.cellmarkers1.norm_rmvfilt[,c(3:12,22:25,28:31)]), center = TRUE, scale = TRUE)
z_dgeFPKM.cellmarkers1.norm_rmvfilt_rmvsamp <- t(z_TdgeFPKM.cellmarkers1.norm_rmvfilt_rmvsamp)

head(z_dgeFPKM.cellmarkers1.norm_rmvfilt_rmvsamp)
dim(z_dgeFPKM.cellmarkers1.norm_rmvfilt_rmvsamp)
#filtered samples
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_dgeFPKM.cellmarkers1.norm_rmvfilt_rmvsamp,
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12,
         filename = "pheatmap_z_dgeFPKM.cellmarkers1.norm_rmvfilt_rmvsamp.png")
dev.off()

#Filtered samples but all genes
pheatmap(z_dgeFPKM.cellmarkers1.norm[,c(3:12,22:25,28:31)],
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12,
         filename = "pheatmap_z_dgeFPKM.cellmarkers1.norm_rmvsamp.png")
dev.off()
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#created after basila gene selection  bouschetmarkerslist.sorted.txt
bouschetmarkerslist.sorted <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/bouschetmarkerslist.sorted.txt", header = F)
head(bouschetmarkerslist.sorted)
colnames(bouschetmarkerslist.sorted) <- c("Gene", "markers")
dim(bouschetmarkerslist.sorted) #print 29 genes
ens_gene_names_chrpos_dedup_M20 <- read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt")
head(ens_gene_names_chrpos_dedup_M20)
ens_gene_names_chrpos_dedup_M20 <- ens_gene_names_chrpos_dedup_M20[,c(3,4,5,1,2)]
head(ens_gene_names_chrpos_dedup_M20)
colnames(ens_gene_names_chrpos_dedup_M20) <- c("chr", "start", "end", "id", "Gene")

bouschetmarkers = merge(bouschetmarkerslist.sorted, ens_gene_names_chrpos_dedup_M20, by="Gene", all.x=FALSE)
head(bouschetmarkers)
dim(bouschetmarkers) #print all 29 genes assigned
bouschetmarkersre <- bouschetmarkers[order(bouschetmarkers$markers),]
head(bouschetmarkersre)
dim(bouschetmarkersre)
#bouschetmarkers.cords <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/bouschetmarkers.symbol.ens.chr.txt", header = F)
#head(bouschetmarkers.cords)
#dim(bouschetmarkers.cords)
#colnames(bouschetmarkers.cords) <- c("chr", "start", "end", "id", "Gene")

dgeFPKM.bouschetmarker = merge(dgeFPKM, bouschetmarkersre, by="id", all.x=FALSE)
head(dgeFPKM.bouschetmarker)
dim(dgeFPKM.bouschetmarker) #Print all 29 genes overlapped and assigned
dgeFPKM.bouschetmarkers <- dgeFPKM.bouschetmarker[order(dgeFPKM.bouschetmarker$markers),]
dim(dgeFPKM.bouschetmarkers)
head(dgeFPKM.bouschetmarkers)
dgeFPKM.bouschetmarkersre <- dgeFPKM.bouschetmarkers[,c(33, 2:32)]
head(dgeFPKM.bouschetmarkersre,2)
rownames(dgeFPKM.bouschetmarkersre) <- dgeFPKM.bouschetmarkersre[,1]
dgeFPKM.bouschetmarkers1 <- dgeFPKM.bouschetmarkersre[,c(2:32)]
head(dgeFPKM.bouschetmarkers1,2)
dgeFPKM.bouschetmarkers1 <- as.matrix(dgeFPKM.bouschetmarkers1)

#Loading control
loadingcontrol.cords <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/loadingcontrol.symbol.ens.chr.txt", header = F)
head(loadingcontrol.cords)
colnames(loadingcontrol.cords) <- c("chr", "start", "end", "id", "Gene")
dgeFPKM.loadingcontrol = merge(dgeFPKM, loadingcontrol.cords, by="id", all.x=FALSE)
head(dgeFPKM.loadingcontrol)
dim(dgeFPKM.loadingcontrol)
dgeFPKM.loadingcontrol <- dgeFPKM.loadingcontrol[,c(33, 2:32)]
head(dgeFPKM.loadingcontrol,2)
rownames(dgeFPKM.loadingcontrol) <- dgeFPKM.loadingcontrol$Gene
dgeFPKM.loadingcontrol1 <- dgeFPKM.loadingcontrol[,c(2:32)]
head(dgeFPKM.loadingcontrol1,2)
dgeFPKM.loadingcontrol1 <- as.matrix(dgeFPKM.loadingcontrol1)
dgeFPKM.loadingcontrol1.Actb <- dgeFPKM.loadingcontrol1[1,]
#Normalize with loading control
dim(dgeFPKM.bouschetmarkers1)

dgeFPKM.bouschetmarkers1.norm <- cbind.data.frame((dgeFPKM.bouschetmarkers1[,1] / dgeFPKM.loadingcontrol1.Actb[1]),
                                                  (dgeFPKM.bouschetmarkers1[,2] / dgeFPKM.loadingcontrol1.Actb[2]),
                                                  (dgeFPKM.bouschetmarkers1[,3] / dgeFPKM.loadingcontrol1.Actb[3]),
                                                  (dgeFPKM.bouschetmarkers1[,4] / dgeFPKM.loadingcontrol1.Actb[4]),
                                                  (dgeFPKM.bouschetmarkers1[,5] / dgeFPKM.loadingcontrol1.Actb[5]),
                                                  (dgeFPKM.bouschetmarkers1[,6] / dgeFPKM.loadingcontrol1.Actb[6]),
                                                  (dgeFPKM.bouschetmarkers1[,7] / dgeFPKM.loadingcontrol1.Actb[7]),
                                                  (dgeFPKM.bouschetmarkers1[,8] / dgeFPKM.loadingcontrol1.Actb[8]),
                                                  (dgeFPKM.bouschetmarkers1[,9] / dgeFPKM.loadingcontrol1.Actb[9]),
                                                  (dgeFPKM.bouschetmarkers1[,10] / dgeFPKM.loadingcontrol1.Actb[10]),
                                                  (dgeFPKM.bouschetmarkers1[,11] / dgeFPKM.loadingcontrol1.Actb[11]),
                                                  (dgeFPKM.bouschetmarkers1[,12] / dgeFPKM.loadingcontrol1.Actb[12]),
                                                  (dgeFPKM.bouschetmarkers1[,13] / dgeFPKM.loadingcontrol1.Actb[13]),
                                                  (dgeFPKM.bouschetmarkers1[,14] / dgeFPKM.loadingcontrol1.Actb[14]),
                                                  (dgeFPKM.bouschetmarkers1[,15] / dgeFPKM.loadingcontrol1.Actb[15]),
                                                  (dgeFPKM.bouschetmarkers1[,16] / dgeFPKM.loadingcontrol1.Actb[16]),
                                                  (dgeFPKM.bouschetmarkers1[,17] / dgeFPKM.loadingcontrol1.Actb[17]),
                                                  (dgeFPKM.bouschetmarkers1[,18] / dgeFPKM.loadingcontrol1.Actb[18]),
                                                  (dgeFPKM.bouschetmarkers1[,19] / dgeFPKM.loadingcontrol1.Actb[19]),
                                                  (dgeFPKM.bouschetmarkers1[,20] / dgeFPKM.loadingcontrol1.Actb[20]),
                                                  (dgeFPKM.bouschetmarkers1[,21] / dgeFPKM.loadingcontrol1.Actb[21]),
                                                  (dgeFPKM.bouschetmarkers1[,22] / dgeFPKM.loadingcontrol1.Actb[22]),
                                                  (dgeFPKM.bouschetmarkers1[,23] / dgeFPKM.loadingcontrol1.Actb[23]),
                                                  (dgeFPKM.bouschetmarkers1[,24] / dgeFPKM.loadingcontrol1.Actb[24]),
                                                  (dgeFPKM.bouschetmarkers1[,25] / dgeFPKM.loadingcontrol1.Actb[25]),
                                                  (dgeFPKM.bouschetmarkers1[,26] / dgeFPKM.loadingcontrol1.Actb[26]),
                                                  (dgeFPKM.bouschetmarkers1[,27] / dgeFPKM.loadingcontrol1.Actb[27]),
                                                  (dgeFPKM.bouschetmarkers1[,28] / dgeFPKM.loadingcontrol1.Actb[28]),
                                                  (dgeFPKM.bouschetmarkers1[,29] / dgeFPKM.loadingcontrol1.Actb[29]),
                                                  (dgeFPKM.bouschetmarkers1[,30] / dgeFPKM.loadingcontrol1.Actb[30]),
                                                  (dgeFPKM.bouschetmarkers1[,31] / dgeFPKM.loadingcontrol1.Actb[31]))
colnames(dgeFPKM.bouschetmarkers1.norm) <- c("ourT0WTE14","ourT0ZFP57KOE14","ourT0WTJB1","ourT0ZFP57KOJB1","SRR1409932","SRR1409933","SRR1409934","SRR1409935","SRR1409938","SRR1409939","SRR1409940","SRR1409941","SRR1409944","SRR1409945","SRR1409946","SRR1409950","SRR1409951","SRR1409952","SRR1409953","SRR1409954","SRR1409955","ourT12WT1", "ourT12WT2", "ourT12ZFP57KO1", "ourT12ZFP57KO2","SRR5665889","SRR5665890","SRR8329324","SRR8329325","SRR8329327","SRR8329328")

head(dgeFPKM.bouschetmarkers1.norm)
boxplot(log10(dgeFPKM.bouschetmarkers1.norm), ylim=c(0,10), col = c("blue","blue","blue","blue","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","red","red","red","red","grey","grey","orange","orange","orange","orange"))
#save as dgeFPKM.bouschetmarkers1.norm.sep.svg
library(EDASeq)
plotPCA(as.matrix(dgeFPKM.bouschetmarkers1.norm))
plotPCA(as.matrix(dgeFPKM.bouschetmarkers1.norm), labels=F, col =  c("cyan","red","cyan","red","blue","blue","darkgreen","darkgreen","blue","blue","darkgreen","darkgreen","darkgreen","darkgreen","darkgreen","seagreen","seagreen","seagreen","seagreen","seagreen","seagreen","green","green","magenta","magenta","darkgreen","darkgreen","navy","navy","red","red"))
ggsave("PCA_tdgeFPKM.bouschetmarkers1.norm_scaleT.sep.svg", width=17*1.25, height=12*1.25, units="cm", dpi=96) #scale =T #library(ggplot2)

write.table(dgeFPKM.bouschetmarkers1.norm, "dgeFPKM.bouschetmarkers1.norm.sep.txt", row.names = T, quote = F, append=F)
#Dendrogram cluster
dim(dgeFPKM.bouschetmarkers1.norm)  #29 31
# Compute distances and hierarchical clustering
dgeFPKM.bouschetmarkers1.norm1 <- t(dgeFPKM.bouschetmarkers1.norm)
dgeFPKM.bouschetmarkers1.normscale <- scale(dgeFPKM.bouschetmarkers1.norm1)
dddgeFPKM.bouschetmarkers1.normscale <- dist(as.matrix(dgeFPKM.bouschetmarkers1.normscale), method = "euclidean")
dddgeFPKM.bouschetmarkers1.normscalehc <- hclust(dddgeFPKM.bouschetmarkers1.normscale, method = "ward.D")
plot(dddgeFPKM.bouschetmarkers1.normscalehc)
#save as dddgeFPKM.bouschetmarkers1.normscalehc.sep.svg
head(dgeFPKM.bouschetmarkers1.norm,1)                                                             
colnames(dgeFPKM.bouschetmarkers1.norm) <- colnames(dgeFPKM.bouschetmarkers1)
head(dgeFPKM.bouschetmarkers1.norm, 2)                                                             
z_TdgeFPKM.bouschetmarkers1.norm= scale(t(dgeFPKM.bouschetmarkers1.norm), center = TRUE, scale = TRUE)
z_dgeFPKM.bouschetmarkers1.norm <- t(z_TdgeFPKM.bouschetmarkers1.norm)
colfunc <- colorRampPalette(c("navy","white","firebrick3"))

dim(z_dgeFPKM.bouschetmarkers1.norm)  #29 31
library(pheatmap)
library(RColorBrewer)
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_dgeFPKM.bouschetmarkers1.norm[,c(1:31)],
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12, 
         filename = "pheatmap_z_dgeFPKM.bouschetmarkers1.norm.sep.png")

dev.off()
#Dendrogram cluster
dim(dgeFPKM.bouschetmarkers1.norm)
head(dgeFPKM.bouschetmarkers1.norm)
#Filter samples and filter genes
#fgrep -f ./../../2019/time12/time0-tim12/removegenesbouschetmarkers.txt dgeFPKM.bouschetmarkers1.norm.sep.txt  -v > dgeFPKM.bouschetmarkers1.norm_rmvfilt.txt
# Compute distances and hierarchical clustering
dgeFPKM.bouschetmarkers1.norm_rmvfilt <- read.table("dgeFPKM.bouschetmarkers1.norm_rmvfilt.txt", header = T, row.names = 1)
dim(dgeFPKM.bouschetmarkers1.norm_rmvfilt)
head(dgeFPKM.bouschetmarkers1.norm_rmvfilt)
#Remove samples E14 cell types
TdgeFPKM.bouschetmarkers1.norm_rmvfilt_rmvsamp <- t(dgeFPKM.bouschetmarkers1.norm_rmvfilt[,c(3:12,22:25,28:31)])
dgeFPKM.bouschetmarkers1.norm_rmvfilt_rmvsamp_scale <- scale(TdgeFPKM.bouschetmarkers1.norm_rmvfilt_rmvsamp)

dddgeFPKM.bouschetmarkers1.norm_rmvfilt_rmvsamp_scale <- dist(as.matrix(dgeFPKM.bouschetmarkers1.norm_rmvfilt_rmvsamp_scale), method = "euclidean")
dddgeFPKM.bouschetmarkers1.norm_rmvfilt_rmvsamp_scalehc <- hclust(dddgeFPKM.bouschetmarkers1.norm_rmvfilt_rmvsamp_scale, method = "ward.D")
plot(dddgeFPKM.bouschetmarkers1.norm_rmvfilt_rmvsamp_scalehc)
#save as dddgeFPKM.bouschetmarkers1.norm_rmvfilt_rmvsamp_scalehc.svg


#z-score convert
z_TdgeFPKM.bouschetmarkers1.norm_rmvfilt_rmvsamp= scale(t(dgeFPKM.bouschetmarkers1.norm_rmvfilt[,c(3:12,22:25,28:31)]), center = TRUE, scale = TRUE)
z_dgeFPKM.bouschetmarkers1.norm_rmvfilt_rmvsamp <- t(z_TdgeFPKM.bouschetmarkers1.norm_rmvfilt_rmvsamp)

head(z_dgeFPKM.bouschetmarkers1.norm_rmvfilt_rmvsamp)
dim(z_dgeFPKM.bouschetmarkers1.norm_rmvfilt_rmvsamp)
#filtered samples
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_dgeFPKM.bouschetmarkers1.norm_rmvfilt_rmvsamp,
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12,
         filename = "pheatmap_z_dgeFPKM.bouschetmarkers1.norm_rmvfilt_rmvsamp.png")
dev.off()

#Filtered samples but all genes
pheatmap(z_dgeFPKM.bouschetmarkers1.norm,
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12,
         filename = "pheatmap_z_dgeFPKM.bouschetmarkers1.norm_rmvsamp.png")
dev.off()
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#Batch Effect detected#
#Remove batch effect using Combat#
sample_info <- read.table("sampleinfoadd.txt", header = TRUE)
head(sample_info)
sample_info[,1]
rownames(sample_info)=sample_info[,1]
sample_info = sample_info[,-1]
head(sample_info)

phenobatch <- read.table("phenoadd.txt", header = TRUE)
rownames(phenobatch)
phenobatch[,1]
rownames(phenobatch)=phenobatch[,1]
rownames(phenobatch)
colnames(phenobatch)
phenobatch = phenobatch[,-1]
head(phenobatch)
colnames(phenobatch)
dim(phenobatch)
Batch = data.frame(phenobatch$Batch)
head(Batch)
rownames(Batch) <- colnames(dgeFPKM[,1:31])
head(Batch)
Run = data.frame(phenobatch$Run)
rownames(Run) <- colnames(dgeFPKM[,1:31])
head(Run)
Modcombat = model.matrix(~1, data=Batch)

#BiocManager::install("sva")
library(sva)
head(dgeFPKM)
dgeFPKMfilt <- dgeFPKM[,1:31]
dgeFPKMfilt[dgeFPKMfilt==0] <- NA
head(dgeFPKMfilt)
dgeFPKMfilt1 <-dgeFPKMfilt[complete.cases(dgeFPKMfilt),]
dim(dgeFPKMfilt1)
head(dgeFPKMfilt1)

head(dgeFPKMfilt1,1)
write.table(dgeFPKMfilt1, "dgeFPKMfilt1.txt", quote=F, append = F, row.names = T)
dgeFPKMcor = ComBat(dat=as.matrix(dgeFPKMfilt1), batch=as.numeric(phenobatch$Batch), mod=Modcombat, par.prior=TRUE, prior.plots=FALSE)
head(dgeFPKMcor)
dim(dgeFPKMcor)
dgeFPKMcor <- data.frame(dgeFPKMcor)
dgeFPKMcor["id"] <- rownames(dgeFPKMcor)
dgeFPKMcor.bouschetmarker = merge(dgeFPKMcor, bouschetmarkersre, by="id", all.x=FALSE)
head(dgeFPKMcor.bouschetmarker)
dim(dgeFPKMcor.bouschetmarker) #Print all 16 genes overlapped and assigned
dgeFPKMcor.bouschetmarkers <- dgeFPKMcor.bouschetmarker[order(dgeFPKMcor.bouschetmarker$markers),]
dim(dgeFPKMcor.bouschetmarkers)
head(dgeFPKMcor.bouschetmarkers)
dgeFPKMcor.bouschetmarkersre <- dgeFPKMcor.bouschetmarkers[,c(33, 2:32)]
head(dgeFPKMcor.bouschetmarkersre,2)
rownames(dgeFPKMcor.bouschetmarkersre) <- dgeFPKMcor.bouschetmarkersre[,1]
dgeFPKMcor.bouschetmarkers1 <- dgeFPKMcor.bouschetmarkersre[,c(2:32)]
head(dgeFPKMcor.bouschetmarkers1,2)
dgeFPKMcor.bouschetmarkers1 <- as.matrix(dgeFPKMcor.bouschetmarkers1)

#Loading control
loadingcontrol.cords <- read.table("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/time0-tim12/loadingcontrol.symbol.ens.chr.txt", header = F)
head(loadingcontrol.cords)
colnames(loadingcontrol.cords) <- c("chr", "start", "end", "id", "Gene")
dgeFPKMcor.loadingcontrol = merge(dgeFPKMcor, loadingcontrol.cords, by="id", all.x=FALSE)
head(dgeFPKMcor.loadingcontrol)
dim(dgeFPKMcor.loadingcontrol)
dgeFPKMcor.loadingcontrol <- dgeFPKMcor.loadingcontrol[,c(33, 2:32)]
head(dgeFPKMcor.loadingcontrol,2)
rownames(dgeFPKMcor.loadingcontrol) <- dgeFPKMcor.loadingcontrol$Gene
dgeFPKMcor.loadingcontrol1 <- dgeFPKMcor.loadingcontrol[,c(2:32)]
head(dgeFPKMcor.loadingcontrol1,2)
dgeFPKMcor.loadingcontrol1 <- as.matrix(dgeFPKMcor.loadingcontrol1)
dgeFPKMcor.loadingcontrol1.Actb <- dgeFPKMcor.loadingcontrol1[1,]
#Normalize with loading control
dim(dgeFPKMcor.bouschetmarkers1)

dgeFPKMcor.bouschetmarkers1.norm <- cbind.data.frame((dgeFPKMcor.bouschetmarkers1[,1] / dgeFPKMcor.loadingcontrol1.Actb[1]),
                                                     (dgeFPKMcor.bouschetmarkers1[,2] / dgeFPKMcor.loadingcontrol1.Actb[2]),
                                                     (dgeFPKMcor.bouschetmarkers1[,3] / dgeFPKMcor.loadingcontrol1.Actb[3]),
                                                     (dgeFPKMcor.bouschetmarkers1[,4] / dgeFPKMcor.loadingcontrol1.Actb[4]),
                                                     (dgeFPKMcor.bouschetmarkers1[,5] / dgeFPKMcor.loadingcontrol1.Actb[5]),
                                                     (dgeFPKMcor.bouschetmarkers1[,6] / dgeFPKMcor.loadingcontrol1.Actb[6]),
                                                     (dgeFPKMcor.bouschetmarkers1[,7] / dgeFPKMcor.loadingcontrol1.Actb[7]),
                                                     (dgeFPKMcor.bouschetmarkers1[,8] / dgeFPKMcor.loadingcontrol1.Actb[8]),
                                                     (dgeFPKMcor.bouschetmarkers1[,9] / dgeFPKMcor.loadingcontrol1.Actb[9]),
                                                     (dgeFPKMcor.bouschetmarkers1[,10] / dgeFPKMcor.loadingcontrol1.Actb[10]),
                                                     (dgeFPKMcor.bouschetmarkers1[,11] / dgeFPKMcor.loadingcontrol1.Actb[11]),
                                                     (dgeFPKMcor.bouschetmarkers1[,12] / dgeFPKMcor.loadingcontrol1.Actb[12]),
                                                     (dgeFPKMcor.bouschetmarkers1[,13] / dgeFPKMcor.loadingcontrol1.Actb[13]),
                                                     (dgeFPKMcor.bouschetmarkers1[,14] / dgeFPKMcor.loadingcontrol1.Actb[14]),
                                                     (dgeFPKMcor.bouschetmarkers1[,15] / dgeFPKMcor.loadingcontrol1.Actb[15]),
                                                     (dgeFPKMcor.bouschetmarkers1[,16] / dgeFPKMcor.loadingcontrol1.Actb[16]),
                                                     (dgeFPKMcor.bouschetmarkers1[,17] / dgeFPKMcor.loadingcontrol1.Actb[17]),
                                                     (dgeFPKMcor.bouschetmarkers1[,18] / dgeFPKMcor.loadingcontrol1.Actb[18]),
                                                     (dgeFPKMcor.bouschetmarkers1[,19] / dgeFPKMcor.loadingcontrol1.Actb[19]),
                                                     (dgeFPKMcor.bouschetmarkers1[,20] / dgeFPKMcor.loadingcontrol1.Actb[20]),
                                                     (dgeFPKMcor.bouschetmarkers1[,21] / dgeFPKMcor.loadingcontrol1.Actb[21]),
                                                     (dgeFPKMcor.bouschetmarkers1[,22] / dgeFPKMcor.loadingcontrol1.Actb[22]),
                                                     (dgeFPKMcor.bouschetmarkers1[,23] / dgeFPKMcor.loadingcontrol1.Actb[23]),
                                                     (dgeFPKMcor.bouschetmarkers1[,24] / dgeFPKMcor.loadingcontrol1.Actb[24]),
                                                     (dgeFPKMcor.bouschetmarkers1[,25] / dgeFPKMcor.loadingcontrol1.Actb[25]),
                                                     (dgeFPKMcor.bouschetmarkers1[,26] / dgeFPKMcor.loadingcontrol1.Actb[26]),
                                                     (dgeFPKMcor.bouschetmarkers1[,27] / dgeFPKMcor.loadingcontrol1.Actb[27]),
                                                     (dgeFPKMcor.bouschetmarkers1[,28] / dgeFPKMcor.loadingcontrol1.Actb[28]),
                                                     (dgeFPKMcor.bouschetmarkers1[,29] / dgeFPKMcor.loadingcontrol1.Actb[29]),
                                                     (dgeFPKMcor.bouschetmarkers1[,30] / dgeFPKMcor.loadingcontrol1.Actb[30]),
                                                     (dgeFPKMcor.bouschetmarkers1[,31] / dgeFPKMcor.loadingcontrol1.Actb[31]))
colnames(dgeFPKMcor.bouschetmarkers1.norm) <- c("ourT0WTE14","ourT0ZFP57KOE14","ourT0WTJB1","ourT0ZFP57KOJB1","SRR1409932","SRR1409933","SRR1409934","SRR1409935","SRR1409938","SRR1409939","SRR1409940","SRR1409941","SRR1409944","SRR1409945","SRR1409946","SRR1409950","SRR1409951","SRR1409952","SRR1409953","SRR1409954","SRR1409955","ourT12WT1", "ourT12WT2", "ourT12ZFP57KO1", "ourT12ZFP57KO2","SRR5665889","SRR5665890","SRR8329324","SRR8329325","SRR8329327","SRR8329328")

head(dgeFPKMcor.bouschetmarkers1.norm)
write.table(dgeFPKMcor.bouschetmarkers1.norm,"dgeFPKMcor.bouschetmarkers1.norm.sep.txt", quote = F, append = F)
#Filter samples and filter genes
#fgrep -f ./../../2019/time12/time0-tim12/removegenesbouschetmarkers.txt dgeFPKMcor.bouschetmarkers1.norm.sep.txt  -v > dgeFPKMcor.bouschetmarkers1.norm_rmvfilt.txt

# Compute distances and hierarchical clustering
dgeFPKMcor.bouschetmarkers1.norm_rmvfilt <- read.table("dgeFPKMcor.bouschetmarkers1.norm_rmvfilt.txt", header = T, row.names = 1)
dim(dgeFPKMcor.bouschetmarkers1.norm_rmvfilt)
head(dgeFPKMcor.bouschetmarkers1.norm_rmvfilt)
##Basilia Suggestion:Hybrid Type1 category
#Remove samples E14 cell types, ferguson E14 cell types
TdgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type1 <- t(dgeFPKMcor.bouschetmarkers1.norm_rmvfilt[,c(3:12,22:25)])
dgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type1_scale <- scale(TdgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type1)

dddgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type1_scale <- dist(as.matrix(dgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type1_scale), method = "euclidean")
dddgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type1_scalehc <- hclust(dddgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type1_scale, method = "ward.D")
plot(dddgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type1_scalehc)
#save as dddgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type1_scalehc.svg


#z-score convert
z_TdgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type1= scale(t(dgeFPKMcor.bouschetmarkers1.norm_rmvfilt[,c(3:12,22:25)]), center = TRUE, scale = TRUE)
z_dgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type1 <- t(z_TdgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type1)

head(z_dgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type1)
dim(z_dgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type1)
#filtered samples
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_dgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type1,
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12,
         filename = "pheatmap_z_dgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type1.png")
dev.off()

##Basilia Suggestion: type2 category
#keep samples E14 cell types, ferguson E14 cell types and in vitro E14 bouschet
TdgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type2 <- t(dgeFPKMcor.bouschetmarkers1.norm_rmvfilt[,c(1:15,22:25,28:31)])
dgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type2_scale <- scale(TdgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type2)

dddgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type2_scale <- dist(as.matrix(dgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type2_scale), method = "euclidean")
dddgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type2_scalehc <- hclust(dddgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type2_scale, method = "ward.D")
plot(dddgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type2_scalehc)
#save as dddgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type2_scalehc.svg


#z-score convert
z_TdgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type2= scale(t(dgeFPKMcor.bouschetmarkers1.norm_rmvfilt[,c(1:15,22:25,28:31)]), center = TRUE, scale = TRUE)
z_dgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type2 <- t(z_TdgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type2)

head(z_dgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type2)
dim(z_dgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type2)
#filtered samples
breaksList1 = seq(-1, 1, by = 0.01)
pheatmap(z_dgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type2,
         color = colorRampPalette(c("navy","white","firebrick3"))(length(breaksList1)),
         breaks = breaksList1,
         fontsize = 8,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = NA,
         cluster_cols=T, 
         cluster_rows=T,
         cutree_cols = 2,
         cellwidth = 20, 
         cellheight = 12,
         filename = "pheatmap_z_dgeFPKMcor.bouschetmarkers1.norm_rmvfilt_type2.png")
dev.off()

################################## Deseq2 check

###################################################### DESeq2 #######################################################
dim(bhybcountdata1)
head(bhybcountdata1)
library(DESeq2)
bhybcountdata1 = as.matrix(bhybcountdata1)
head(bhybcountdata1)
coldatarisot0 <- read.table("coldatarisot0.txt" , header = TRUE, stringsAsFactors = FALSE)
rownames(coldatarisot0)
coldatarisot0[,1]
rownames(coldatarisot0)=coldatarisot0[,1]
rownames(coldatarisot0)
colnames(coldatarisot0)
coldatarisot0 = coldatarisot0[,-1]
coldatarisot0 <- data.frame(coldatarisot0)
head(coldatarisot0)
rownames(coldatarisot0) <- colnames(bhybcountdata1)
coldatarisot0 <- coldatarisot0[,c("condition","replicate")]
coldatarisot0$condition <- factor(coldatarisot0$condition)
coldatarisot0$replicate <- factor(coldatarisot0$replicate)
all(rownames(coldatarisot0) == colnames(bhybcountdata1)) #should print TRUE
ddsrisot0 <- DESeqDataSetFromMatrix(countData =bhybcountdata1, colData = coldatarisot0, design = ~ condition)


ddsrisot0
featureData <- data.frame(gene=rownames(bhybcountdata1))
keep <- rowSums(counts(ddsrisot0)) >= 9
ddsrisot0 <- ddsrisot0[keep,]

#View filtered count matrix: View(counts(ddsrisot0))
#Normalization is the part of DESeq command: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#Normalized separately:
ddsrisot0Norm <- estimateSizeFactors(ddsrisot0)
sizeFactors(ddsrisot0Norm)
#Export normalized counts
ddsrisot0Normcounts <- counts(ddsrisot0Norm, normalized=TRUE)
head(ddsrisot0Normcounts)
write.table(ddsrisot0Normcounts, "ddsrisot0Normcounts.txt", sep="\t", quote=F, col.names=NA)

ddsrisot0Normcounts <- data.frame(ddsrisot0Normcounts)
#Chromosome positions
chr.pos =  read.table("/home/ankitv/ref_av/gencodes/gencode_M20/prep/ens_gene_names_chrpos_dedup_M20.txt",header=FALSE)
head(chr.pos)
colnames(chr.pos) <- c("ensid", "Genes", "chr", "start", "end")
head(chr.pos)
head(ddsrisot0Normcounts)
ddsrisot0Normcounts["ensid"] <- rownames(ddsrisot0Normcounts)
ddsrisot0Normcounts_chr.pos = merge(ddsrisot0Normcounts, chr.pos, by="ensid", all.x=TRUE)
head(ddsrisot0Normcounts_chr.pos)
dim(ddsrisot0Normcounts_chr.pos)
ddsrisot0Normcounts_chr.pos <- ddsrisot0Normcounts_chr.pos[,c(5:7,4,1:3)]
write.table(ddsrisot0Normcounts_chr.pos, "ddsrisot0Normcounts_chr.pos_gene.txt", sep="\t", quote = FALSE, append = FALSE)
