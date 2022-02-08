#~/tools_av/STAR-2.5.4a/bin/Linux_x86_64/STAR --runThreadN 10 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles mm10.fa --sjdbGTFfile gencode.vM20.chr_patch_hapl_scaff.annotation.gtf --sjdbOverhang 124
#gencode.vM20.chr_patch_hapl_scaff.annotation.gtf contains format "chr" and mm10.fa also contains >chr*.fa : both  have chr
cd JB1_WT
cd rep1
~/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode alignReads --runThreadN 10 --genomeDir /home/ankitv/ref_av/gencodes/gencode_M20/ --sjdbGTFfile ~/ref_av/gencodes/gencode_M20/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf --sjdbOverhang 124 --readFilesIn /home/ankitv/2019/rna-seq/time12/JB1_WT/rep1/JB1_WT_Rep1_R1_paired.fastq.gz /home/ankitv/2019/rna-seq/time12/JB1_WT/rep1/JB1_WT_Rep1_R2_paired.fastq.gz --outFileNamePrefix JB1_WT_Rep1_ --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outBAMsortingThreadN 12 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMattributes All --outSAMstrandField intronMotif --readFilesCommand gunzip -c
samtools index *.sortedByCoord.out.bam
/home/ankitv/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --runThreadN 10 --inputBAMfile JB1_WT_Rep1_Aligned.sortedByCoord.out.bam --outWigType bedGraph --outWigNorm RPM --outWigStrand Stranded --outFileNamePrefix JB1_WT_Rep1_time12_All_
grep "NT" JB1_WT_Rep1_time12_All_Signal.Unique.str1.out.bg -v > JB1_WT_Rep1_time12_All_Signal.Unique.str1.out-NT.bg
sort -k1,1 -k2,2n JB1_WT_Rep1_time12_All_Signal.Unique.str1.out-NT.bg > JB1_WT_Rep1_time12_All_Signal.Unique.str1.out.-NT.sort.bg
grep "NT" JB1_WT_Rep1_time12_All_Signal.Unique.str2.out.bg -v > JB1_WT_Rep1_time12_All_Signal.Unique.str2.out-NT.bg
sort -k1,1 -k2,2n JB1_WT_Rep1_time12_All_Signal.Unique.str2.out-NT.bg > JB1_WT_Rep1_time12_All_Signal.Unique.str2.out.-NT.sort.bg
~/tools_av/bedGraphToBigWig JB1_WT_Rep1_time12_All_Signal.Unique.str1.out.-NT.sort.bg ~/ref_av/mm10/mm10.chrom.sizes All_T12_JB1_WT_Rep1_str1_mm10.bw
~/tools_av/bedGraphToBigWig JB1_WT_Rep1_time12_All_Signal.Unique.str2.out.-NT.sort.bg ~/ref_av/mm10/mm10.chrom.sizes All_T12_JB1_WT_Rep1_str2_mm10.bw
rm -r JB1_WT_Rep1_time12_All_Signal.Unique.str1.out-NT.bg
rm -r JB1_WT_Rep1_time12_All_Signal.Unique.str2.out-NT.bg
gzip *-NT.sort.bg
cd ..
cd rep2
~/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode alignReads --runThreadN 10 --genomeDir /home/ankitv/ref_av/gencodes/gencode_M20/ --sjdbGTFfile ~/ref_av/gencodes/gencode_M20/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf --sjdbOverhang 124 --readFilesIn /home/ankitv/2019/rna-seq/time12/JB1_WT/rep2/JB1_WT_Rep2_R1_paired.fastq.gz /home/ankitv/2019/rna-seq/time12/JB1_WT/rep2/JB1_WT_Rep2_R2_paired.fastq.gz --outFileNamePrefix JB1_WT_Rep2_ --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outBAMsortingThreadN 12 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMattributes All --outSAMstrandField intronMotif --readFilesCommand gunzip -c
samtools index *.sortedByCoord.out.bam
/home/ankitv/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --runThreadN 12 --inputBAMfile JB1_WT_Rep2_Aligned.sortedByCoord.out.bam --outWigType bedGraph --outWigNorm RPM --outWigStrand Stranded --outFileNamePrefix JB1_WT_Rep2_time12_All_
grep "NT" JB1_WT_Rep2_time12_All_Signal.Unique.str1.out.bg -v > JB1_WT_Rep2_time12_All_Signal.Unique.str1.out-NT.bg
sort -k1,1 -k2,2n JB1_WT_Rep2_time12_All_Signal.Unique.str1.out-NT.bg > JB1_WT_Rep2_time12_All_Signal.Unique.str1.out.-NT.sort.bg
grep "NT" JB1_WT_Rep2_time12_All_Signal.Unique.str2.out.bg -v > JB1_WT_Rep2_time12_All_Signal.Unique.str2.out-NT.bg
sort -k1,1 -k2,2n JB1_WT_Rep2_time12_All_Signal.Unique.str2.out-NT.bg > JB1_WT_Rep2_time12_All_Signal.Unique.str2.out.-NT.sort.bg
~/tools_av/bedGraphToBigWig JB1_WT_Rep2_time12_All_Signal.Unique.str1.out.-NT.sort.bg ~/ref_av/mm10/mm10.chrom.sizes All_T12_JB1_WT_Rep2_str1_mm10.bw
~/tools_av/bedGraphToBigWig JB1_WT_Rep2_time12_All_Signal.Unique.str2.out.-NT.sort.bg ~/ref_av/mm10/mm10.chrom.sizes All_T12_JB1_WT_Rep2_str2_mm10.bw
rm -r JB1_WT_Rep2_time12_All_Signal.Unique.str1.out-NT.bg
rm -r JB1_WT_Rep2_time12_All_Signal.Unique.str2.out-NT.bg
gzip *-NT.sort.bg
cd ..
cd ..
cd JB1_ZFP57_KO
cd rep1
~/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode alignReads --runThreadN 10 --genomeDir /home/ankitv/ref_av/gencodes/gencode_M20/ --sjdbGTFfile ~/ref_av/gencodes/gencode_M20/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf --sjdbOverhang 124 --readFilesIn /home/ankitv/2019/rna-seq/time12/JB1_ZFP57_KO/rep1/JB1_ZFP57_KO_Rep1_R1_paired.fastq.gz /home/ankitv/2019/rna-seq/time12/JB1_ZFP57_KO/rep1/JB1_ZFP57_KO_Rep1_R2_paired.fastq.gz --outFileNamePrefix JB1_ZFP57_KO_Rep1_ --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outBAMsortingThreadN 12 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMattributes All --outSAMstrandField intronMotif --readFilesCommand gunzip -c
samtools index *.sortedByCoord.out.bam
/home/ankitv/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --runThreadN 12 --inputBAMfile JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bam --outWigType bedGraph --outWigNorm RPM --outWigStrand Stranded --outFileNamePrefix JB1_ZFP57_KO_Rep1_time12_All_
grep "NT" JB1_ZFP57_KO_Rep1_time12_All_Signal.Unique.str1.out.bg -v > JB1_ZFP57_KO_Rep1_time12_All_Signal.Unique.str1.out-NT.bg
sort -k1,1 -k2,2n JB1_ZFP57_KO_Rep1_time12_All_Signal.Unique.str1.out-NT.bg > JB1_ZFP57_KO_Rep1_time12_All_Signal.Unique.str1.out.-NT.sort.bg
grep "NT" JB1_ZFP57_KO_Rep1_time12_All_Signal.Unique.str2.out.bg -v > JB1_ZFP57_KO_Rep1_time12_All_Signal.Unique.str2.out-NT.bg
sort -k1,1 -k2,2n JB1_ZFP57_KO_Rep1_time12_All_Signal.Unique.str2.out-NT.bg > JB1_ZFP57_KO_Rep1_time12_All_Signal.Unique.str2.out.-NT.sort.bg
~/tools_av/bedGraphToBigWig JB1_ZFP57_KO_Rep1_time12_All_Signal.Unique.str1.out.-NT.sort.bg ~/ref_av/mm10/mm10.chrom.sizes All_T12_JB1_ZFP57_KO_Rep1_str1_mm10.bw
~/tools_av/bedGraphToBigWig JB1_ZFP57_KO_Rep1_time12_All_Signal.Unique.str2.out.-NT.sort.bg ~/ref_av/mm10/mm10.chrom.sizes All_T12_JB1_ZFP57_KO_Rep1_str2_mm10.bw
rm -r JB1_ZFP57_KO_Rep1_time12_All_Signal.Unique.str1.out-NT.bg
rm -r JB1_ZFP57_KO_Rep1_time12_All_Signal.Unique.str2.out-NT.bg
gzip *-NT.sort.bg
cd ..
cd rep2
~/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode alignReads --runThreadN 10 --genomeDir /home/ankitv/ref_av/gencodes/gencode_M20/ --sjdbGTFfile ~/ref_av/gencodes/gencode_M20/gencode.vM20.chr_patch_hapl_scaff.annotation.gtf --sjdbOverhang 124 --readFilesIn /home/ankitv/2019/rna-seq/time12/JB1_ZFP57_KO/rep2/JB1_ZFP57_KO_Rep2_R1_paired.fastq.gz /home/ankitv/2019/rna-seq/time12/JB1_ZFP57_KO/rep2/JB1_ZFP57_KO_Rep2_R2_paired.fastq.gz --outFileNamePrefix JB1_ZFP57_KO_Rep2_ --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outBAMsortingThreadN 12 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMattributes All --outSAMstrandField intronMotif --readFilesCommand gunzip -c
samtools index *.sortedByCoord.out.bam
/home/ankitv/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --runThreadN 12 --inputBAMfile JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bam --outWigType bedGraph --outWigNorm RPM --outWigStrand Stranded --outFileNamePrefix JB1_ZFP57_KO_Rep2_time12_All_
grep "NT" JB1_ZFP57_KO_Rep2_time12_All_Signal.Unique.str1.out.bg -v > JB1_ZFP57_KO_Rep2_time12_All_Signal.Unique.str1.out-NT.bg
sort -k1,1 -k2,2n JB1_ZFP57_KO_Rep2_time12_All_Signal.Unique.str1.out-NT.bg > JB1_ZFP57_KO_Rep2_time12_All_Signal.Unique.str1.out.-NT.sort.bg
grep "NT" JB1_ZFP57_KO_Rep2_time12_All_Signal.Unique.str2.out.bg -v > JB1_ZFP57_KO_Rep2_time12_All_Signal.Unique.str2.out-NT.bg
sort -k1,1 -k2,2n JB1_ZFP57_KO_Rep2_time12_All_Signal.Unique.str2.out-NT.bg > JB1_ZFP57_KO_Rep2_time12_All_Signal.Unique.str2.out.-NT.sort.bg
~/tools_av/bedGraphToBigWig JB1_ZFP57_KO_Rep2_time12_All_Signal.Unique.str1.out.-NT.sort.bg ~/ref_av/mm10/mm10.chrom.sizes All_T12_JB1_ZFP57_KO_Rep2_str1_mm10.bw
~/tools_av/bedGraphToBigWig JB1_ZFP57_KO_Rep2_time12_All_Signal.Unique.str2.out.-NT.sort.bg ~/ref_av/mm10/mm10.chrom.sizes All_T12_JB1_ZFP57_KO_Rep2_str2_mm10.bw
rm -r JB1_ZFP57_KO_Rep2_time12_All_Signal.Unique.str1.out-NT.bg
rm -r JB1_ZFP57_KO_Rep2_time12_All_Signal.Unique.str2.out-NT.bg
gzip *-NT.sort.bg
cd ..
cd ..

#sorted by readname
samtools sort -n JB1_WT_Rep1_Aligned.sortedByCoord.out.bam -o JB1_WT_Rep1_Aligned.sortedByReadname.out.bam
samtools sort -n JB1_WT_Rep2_Aligned.sortedByCoord.out.bam -o JB1_WT_Rep2_Aligned.sortedByReadname.out.bam
samtools sort -n JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bam -o JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.bam
samtools sort -n JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bam -o JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.bam

#The scaling factor estimated from Bulk data is as follows by samtools flagstat and will be used for all allele specific bigwig conversion
#JB1_WT_Rep1 = 0.03895
#JB1_WT_Rep2 = 0.03872
#JB1_ZFP57KO_Rep1 =0.03863 
#JB1_ZFP57KO_Rep2 =0.03785

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
