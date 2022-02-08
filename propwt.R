
#------------  Proportion test -------------#
#Combined p.vales -2Sumation(log(Pi))
#ChiSquare df = 2*length(p-values), so if 2 p.values 2 *2
#lower.tail =FALSE and lower.tail =TRUE, 1-p
#Prop test and tag chromsome positions
setwd("/media/ankitv/Archivio2/ankit/rna-seq/mouse/2019/time12/allele_specific_mm10/GRCm38_mm10/allele_sp/featurecount/5Mb/propwt")

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
#Filter B6+JF1 > 10 in WT, if both the alleles are 10 in either of the replicates in WT
gddsaSFNormcounts_chr_gene1_prop.filt <- gddsaSFNormcounts_chr_gene1_prop[which(gddsaSFNormcounts_chr_gene1_prop$JB1_WT_alavg.Rep1 > 10 &
                                                                                  gddsaSFNormcounts_chr_gene1_prop$JB1_WT_alavg.Rep2 > 10),]


dim(gddsaSFNormcounts_chr_gene1_prop.filt)#11217    18
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
imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af <- imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher[which(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher$WT_allele_ratio_1 <= 0.33 |
                                                                                                                                                                                                                                      imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher$WT_allele_ratio_1 >= 0.67),]

dim(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af)

imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af2 <- imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af[which(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af$WT_allele_ratio_2 <= 0.33 |
                                                                                                                                                                                                                                          imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af$WT_allele_ratio_2 >= 0.67),]
dim(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af2)
#Remove very low expresed genes
imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af3 <- imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af2[which(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af2$Gene != "B830012L14Rik"),]

dim(imprinted_gddsaSFNormcounts_chr_gene1_prop.filt_gene_WT_fdrsig_minusX_JB1_Zfp57_overlapped_Kap1_mm10_fisher_af3)
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

#Coloring
fgrep -f imprinted_gene_name.txt  Allele_Specific_Supplementary_Table_4.txt -w > Allele_Specific_Supplementary_Table_4_imprinted_genes.txt
fgrep -f imprinted_gene_name.txt  Allele_Specific_Supplementary_Table_4.txt -w -v > Allele_Specific_Supplementary_Table_4_Nonimprinted_genes.txt
#Open one by one and paste in excel sheet
#First paste imprinted color red all genes
#Now paste non-imprinted and take header and put it in first line
#Benefit of doing this it will also take care of NA and mention blank as NA
#Color by alleles, sort by chr and start---Manually checked the sorting if done proterly
#Name columns as required, 
#Also in excel sheet manually filtered NA from adj p-value, leftover is 11371 which is equal to dataframe used for prop test
#data in excel test_roughNA with NA shows either B6+JF1 <10 in WT 
#save as Allele_Specific_Supplementary_Table_4_Arranged.xlsx and send to Acurzio. 





#Prepare Alsp Table all remianing genes after filtered and prop.test 
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
#save as Allele_Specific_ALL_Supplementary_Table_4_Arranged.xlsx and send to Acurzio. She has to sort by chromosme and start positions wise

########################  End of Analysis ##################################à

