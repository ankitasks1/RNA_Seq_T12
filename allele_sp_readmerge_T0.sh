#The scaling factor estimated from Bulk data is as follows by samtools flagstat and will be used for all allele specific bigwig conversion
#		Aligned	JB1_WT	JB1_ZFP57KO
#		Rep1	9094564	12514695
#		Rep2	9890390	13641894
#		Rep3	28383708	38412364
#			
#		SUM	47368662	64568953
#		RPM	0.0211110037	0.0154873194
#Use for bedgraph normalization		Roundoff	0.0211	0.0155


cd JB1_WT/
~/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode alignReads --alignEndsType EndToEnd --runThreadN 12 --genomeDir /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1_GRCm38_50_overlapped/ --sjdbGTFfile /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1_GRCm38_50_overlapped/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschr.gtf --sjdbOverhang 49 --readFilesIn /media/ankitv/Archivio2/ankit/rna-seq/mouse/2018/time0/2015_05_28_JB1_WT_Zfp57--_clone4/DATI_GREZZI/JB1_WT/JB1_WT_trimmed.fastq.gz --outFileNamePrefix JB1_WT_ --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outBAMsortingThreadN 12 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMattributes All --outSAMstrandField intronMotif --readFilesCommand gunzip -c


samtools sort -n JB1_WT_Aligned.sortedByCoord.out.bam -o JB1_WT_Aligned.sortedByReadname.out.bam
~/tools_av/SNPsplit_v0.3.2/SNPsplit --snp_file /home/ankitv/ref_av/gencodes/gencode_M20/jf1v2_Snp-chr.GRCm38.mm10.Snpfile.txt JB1_WT_Aligned.sortedByReadname.out.bam
samtools sort -o JB1_WT_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_WT_Aligned.sortedByReadname.out.genome1.bam 
samtools sort -o JB1_WT_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_WT_Aligned.sortedByReadname.out.genome2.bam



/home/ankitv/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --runThreadN 12 --inputBAMfile JB1_WT_Aligned.sortedByReadname.out.genome.sort.B6.bam --outWigType bedGraph --outWigNorm None --outWigStrand Stranded --outFileNamePrefix JB1_WT_time0_AlSp_B6_
grep "NT" JB1_WT_time0_AlSp_B6_Signal.Unique.str1.out.bg -v > JB1_WT_time0_AlSp_B6_Signal.Unique.str1.out-NT.bg
sort -k1,1 -k2,2n JB1_WT_time0_AlSp_B6_Signal.Unique.str1.out-NT.bg > JB1_WT_time0_AlSp_B6_Signal.Unique.str1.out.-NT.sort.bg
grep "NT" JB1_WT_time0_AlSp_B6_Signal.Unique.str2.out.bg -v > JB1_WT_time0_AlSp_B6_Signal.Unique.str2.out-NT.bg
sort -k1,1 -k2,2n JB1_WT_time0_AlSp_B6_Signal.Unique.str2.out-NT.bg > JB1_WT_time0_AlSp_B6_Signal.Unique.str2.out.-NT.sort.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4*0.0211}' JB1_WT_time0_AlSp_B6_Signal.Unique.str1.out.-NT.sort.bg > JB1_WT_time0_AlSp_B6_Signal.Unique.str1.out.-NT.sort.chr.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4*0.0211}' JB1_WT_time0_AlSp_B6_Signal.Unique.str2.out.-NT.sort.bg > JB1_WT_time0_AlSp_B6_Signal.Unique.str2.out.-NT.sort.chr.bg
~/tools_av/bedGraphToBigWig JB1_WT_time0_AlSp_B6_Signal.Unique.str1.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_B6_T0_JB1_WT_str1.GRCm38.mm10.bw
~/tools_av/bedGraphToBigWig JB1_WT_time0_AlSp_B6_Signal.Unique.str2.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_B6_T0_JB1_WT_str2.GRCm38.mm10.bw
rm -r JB1_WT_time0_AlSp_B6_Signal.Unique.str1.out-NT.bg
rm -r JB1_WT_time0_AlSp_B6_Signal.Unique.str2.out-NT.bg
/home/ankitv/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --runThreadN 12 --inputBAMfile JB1_WT_Aligned.sortedByReadname.out.genome.sort.JF1.bam --outWigType bedGraph --outWigNorm None --outWigStrand Stranded --outFileNamePrefix JB1_WT_time0_AlSp_JF1_
grep "NT" JB1_WT_time0_AlSp_JF1_Signal.Unique.str1.out.bg -v > JB1_WT_time0_AlSp_JF1_Signal.Unique.str1.out-NT.bg
sort -k1,1 -k2,2n JB1_WT_time0_AlSp_JF1_Signal.Unique.str1.out-NT.bg > JB1_WT_time0_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.bg
grep "NT" JB1_WT_time0_AlSp_JF1_Signal.Unique.str2.out.bg -v > JB1_WT_time0_AlSp_JF1_Signal.Unique.str2.out-NT.bg
sort -k1,1 -k2,2n JB1_WT_time0_AlSp_JF1_Signal.Unique.str2.out-NT.bg > JB1_WT_time0_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4*0.0211}' JB1_WT_time0_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.bg > JB1_WT_time0_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.chr.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4*0.0211}' JB1_WT_time0_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.bg > JB1_WT_time0_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.chr.bg
~/tools_av/bedGraphToBigWig JB1_WT_time0_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_JF1_T0_JB1_WT_str1.GRCm38.mm10.bw
~/tools_av/bedGraphToBigWig JB1_WT_time0_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_JF1_T0_JB1_WT_str2.GRCm38.mm10.bw
rm -r JB1_WT_time0_AlSp_JF1_Signal.Unique.str1.out-NT.bg
rm -r JB1_WT_time0_AlSp_JF1_Signal.Unique.str2.out-NT.bg
gzip *.bg

cd ..
cd ..
cd JB1_ZFP57_KO
cd rep1
~/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode alignReads --alignEndsType EndToEnd --runThreadN 12 --genomeDir /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1_GRCm38_50_overlapped/ --sjdbGTFfile /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1_GRCm38_50_overlapped/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschr.gtf --sjdbOverhang 49 --readFilesIn /media/ankitv/Archivio2/ankit/rna-seq/mouse/2018/time0/2015_05_28_JB1_WT_Zfp57--_clone4/DATI_GREZZI/JB1_ZFP57_KO/JB1_ZFP57_KO_trimmed.fastq.gz --outFileNamePrefix JB1_ZFP57_KO_ --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outBAMsortingThreadN 12 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMattributes All --outSAMstrandField intronMotif --readFilesCommand gunzip -c

samtools sort -n JB1_ZFP57_KO_Aligned.sortedByCoord.out.bam -o JB1_ZFP57_KO_Aligned.sortedByReadname.out.bam


~/tools_av/SNPsplit_v0.3.2/SNPsplit --snp_file /home/ankitv/ref_av/gencodes/gencode_M20/jf1v2_Snp-chr.GRCm38.mm10.Snpfile.txt JB1_ZFP57_KO_Aligned.sortedByReadname.out.bam
samtools sort -o JB1_ZFP57_KO_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Aligned.sortedByReadname.out.genome1.bam 
samtools sort -o JB1_ZFP57_KO_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Aligned.sortedByReadname.out.genome2.bam




/home/ankitv/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --runThreadN 12 --inputBAMfile JB1_ZFP57_KO_Aligned.sortedByReadname.out.genome.sort.B6.bam --outWigType bedGraph --outWigNorm None --outWigStrand Stranded --outFileNamePrefix JB1_ZFP57_KO_time0_AlSp_B6_
grep "NT" JB1_ZFP57_KO_time0_AlSp_B6_Signal.Unique.str1.out.bg -v > JB1_ZFP57_KO_time0_AlSp_B6_Signal.Unique.str1.out-NT.bg
sort -k1,1 -k2,2n JB1_ZFP57_KO_time0_AlSp_B6_Signal.Unique.str1.out-NT.bg > JB1_ZFP57_KO_time0_AlSp_B6_Signal.Unique.str1.out.-NT.sort.bg
grep "NT" JB1_ZFP57_KO_time0_AlSp_B6_Signal.Unique.str2.out.bg -v > JB1_ZFP57_KO_time0_AlSp_B6_Signal.Unique.str2.out-NT.bg
sort -k1,1 -k2,2n JB1_ZFP57_KO_time0_AlSp_B6_Signal.Unique.str2.out-NT.bg > JB1_ZFP57_KO_time0_AlSp_B6_Signal.Unique.str2.out.-NT.sort.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4*0.0155}' JB1_ZFP57_KO_time0_AlSp_B6_Signal.Unique.str1.out.-NT.sort.bg > JB1_ZFP57_KO_time0_AlSp_B6_Signal.Unique.str1.out.-NT.sort.chr.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4*0.0155}' JB1_ZFP57_KO_time0_AlSp_B6_Signal.Unique.str2.out.-NT.sort.bg > JB1_ZFP57_KO_time0_AlSp_B6_Signal.Unique.str2.out.-NT.sort.chr.bg
~/tools_av/bedGraphToBigWig JB1_ZFP57_KO_time0_AlSp_B6_Signal.Unique.str1.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_B6_T0_JB1_ZFP57_KO_str1.GRCm38.mm10.bw
~/tools_av/bedGraphToBigWig JB1_ZFP57_KO_time0_AlSp_B6_Signal.Unique.str2.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_B6_T0_JB1_ZFP57_KO_str2.GRCm38.mm10.bw
rm -r JB1_ZFP57_KO_time0_AlSp_B6_Signal.Unique.str1.out-NT.bg
rm -r JB1_ZFP57_KO_time0_AlSp_B6_Signal.Unique.str2.out-NT.bg
/home/ankitv/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --runThreadN 12 --inputBAMfile JB1_ZFP57_KO_Aligned.sortedByReadname.out.genome.sort.JF1.bam --outWigType bedGraph --outWigNorm None --outWigStrand Stranded --outFileNamePrefix JB1_ZFP57_KO_time0_AlSp_JF1_
grep "NT" JB1_ZFP57_KO_time0_AlSp_JF1_Signal.Unique.str1.out.bg -v > JB1_ZFP57_KO_time0_AlSp_JF1_Signal.Unique.str1.out-NT.bg
sort -k1,1 -k2,2n JB1_ZFP57_KO_time0_AlSp_JF1_Signal.Unique.str1.out-NT.bg > JB1_ZFP57_KO_time0_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.bg
grep "NT" JB1_ZFP57_KO_time0_AlSp_JF1_Signal.Unique.str2.out.bg -v > JB1_ZFP57_KO_time0_AlSp_JF1_Signal.Unique.str2.out-NT.bg
sort -k1,1 -k2,2n JB1_ZFP57_KO_time0_AlSp_JF1_Signal.Unique.str2.out-NT.bg > JB1_ZFP57_KO_time0_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4*0.0155}' JB1_ZFP57_KO_time0_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.bg > JB1_ZFP57_KO_time0_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.chr.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4*0.0155}' JB1_ZFP57_KO_time0_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.bg > JB1_ZFP57_KO_time0_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.chr.bg
~/tools_av/bedGraphToBigWig JB1_ZFP57_KO_time0_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_JF1_T0_JB1_ZFP57_KO_str1.GRCm38.mm10.bw
~/tools_av/bedGraphToBigWig JB1_ZFP57_KO_time0_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_JF1_T0_JB1_ZFP57_KO_str2.GRCm38.mm10.bw
rm -r JB1_ZFP57_KO_time0_AlSp_JF1_Signal.Unique.str1.out-NT.bg
rm -r JB1_ZFP57_KO_time0_AlSp_JF1_Signal.Unique.str2.out-NT.bg
gzip *.bg

cd ..
#igf2, SE gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschrlike-.gtf can be used as it has no readlength isssue
##/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-GRCm38-M20-overlapped/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschrlike-.gtf -T 12 -o Alsp_JB1_time0_star-featureCounts_GRCm38.mm10_Igf2.txt JB1_WT_Aligned.sortedByReadname.out.genome.sort.B6.bam  JB1_WT_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Aligned.sortedByReadname.out.genome.sort.JF1.bam




#/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1_GRCm38_50_overlapped/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschr.gtf -T 12 -o Alsp_JB1_time0_star-featureCounts_GRCm38.txt JB1_WT_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_WT_Rep3_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_WT_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_WT_Rep3_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Rep3_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Rep3_Aligned.sortedByReadname.out.genome.sort.JF1.bam
