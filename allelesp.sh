cd JB1_WT
cd rep1
~/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode alignReads --alignEndsType EndToEnd --runThreadN 10 --genomeDir /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-Grcm38_overlapped/ --sjdbGTFfile /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-Grcm38_overlapped/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschr.gtf --sjdbOverhang 124 --readFilesIn /home/ankitv/2019/rna-seq/time12/JB1_WT/rep1/JB1_WT_Rep1_R1_paired.fastq.gz /home/ankitv/2019/rna-seq/time12/JB1_WT/rep1/JB1_WT_Rep1_R2_paired.fastq.gz --outFileNamePrefix JB1_WT_Rep1_ --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outBAMsortingThreadN 12 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMattributes All --outSAMstrandField intronMotif --readFilesCommand gunzip -c
samtools sort -n JB1_WT_Rep1_Aligned.sortedByCoord.out.bam -o JB1_WT_Rep1_Aligned.sortedByReadname.out.bam
~/tools_av/SNPsplit_v0.3.2/SNPsplit --paired --snp_file /home/ankitv/ref_av/gencodes/gencode_M20/jf1v2_Snp-chr.GRCm38.mm10.Snpfile.txt JB1_WT_Rep1_Aligned.sortedByReadname.out.bam
samtools sort -o JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_WT_Rep1_Aligned.sortedByReadname.out.genome1.bam 
samtools sort -o JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_WT_Rep1_Aligned.sortedByReadname.out.genome2.bam
/home/ankitv/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --runThreadN 12 --inputBAMfile JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bam --outWigType bedGraph --outWigNorm RPM --outWigStrand Stranded --outFileNamePrefix JB1_WT_Rep1_time12_AlSp_B6_
grep "NT" JB1_WT_Rep1_time12_AlSp_B6_Signal.Unique.str1.out.bg -v > JB1_WT_Rep1_time12_AlSp_B6_Signal.Unique.str1.out-NT.bg
sort -k1,1 -k2,2n JB1_WT_Rep1_time12_AlSp_B6_Signal.Unique.str1.out-NT.bg > JB1_WT_Rep1_time12_AlSp_B6_Signal.Unique.str1.out.-NT.sort.bg
grep "NT" JB1_WT_Rep1_time12_AlSp_B6_Signal.Unique.str2.out.bg -v > JB1_WT_Rep1_time12_AlSp_B6_Signal.Unique.str2.out-NT.bg
sort -k1,1 -k2,2n JB1_WT_Rep1_time12_AlSp_B6_Signal.Unique.str2.out-NT.bg > JB1_WT_Rep1_time12_AlSp_B6_Signal.Unique.str2.out.-NT.sort.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' JB1_WT_Rep1_time12_AlSp_B6_Signal.Unique.str1.out.-NT.sort.bg > JB1_WT_Rep1_time12_AlSp_B6_Signal.Unique.str1.out.-NT.sort.chr.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' JB1_WT_Rep1_time12_AlSp_B6_Signal.Unique.str2.out.-NT.sort.bg > JB1_WT_Rep1_time12_AlSp_B6_Signal.Unique.str2.out.-NT.sort.chr.bg
~/tools_av/bedGraphToBigWig JB1_WT_Rep1_time12_AlSp_B6_Signal.Unique.str1.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_B6_T12_JB1_WT_Rep1_str1_GRCm38.mm10.bw
~/tools_av/bedGraphToBigWig JB1_WT_Rep1_time12_AlSp_B6_Signal.Unique.str2.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_B6_T12_JB1_WT_Rep1_str2_GRCm38.mm10.bw
rm -r JB1_WT_Rep1_time12_AlSp_B6_Signal.Unique.str1.out-NT.bg
rm -r JB1_WT_Rep1_time12_AlSp_B6_Signal.Unique.str2.out-NT.bg
/home/ankitv/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --runThreadN 12 --inputBAMfile JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bam --outWigType bedGraph --outWigNorm RPM --outWigStrand Stranded --outFileNamePrefix JB1_WT_Rep1_time12_AlSp_JF1_
grep "NT" JB1_WT_Rep1_time12_AlSp_JF1_Signal.Unique.str1.out.bg -v > JB1_WT_Rep1_time12_AlSp_JF1_Signal.Unique.str1.out-NT.bg
sort -k1,1 -k2,2n JB1_WT_Rep1_time12_AlSp_JF1_Signal.Unique.str1.out-NT.bg > JB1_WT_Rep1_time12_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.bg
grep "NT" JB1_WT_Rep1_time12_AlSp_JF1_Signal.Unique.str2.out.bg -v > JB1_WT_Rep1_time12_AlSp_JF1_Signal.Unique.str2.out-NT.bg
sort -k1,1 -k2,2n JB1_WT_Rep1_time12_AlSp_JF1_Signal.Unique.str2.out-NT.bg > JB1_WT_Rep1_time12_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' JB1_WT_Rep1_time12_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.bg > JB1_WT_Rep1_time12_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.chr.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' JB1_WT_Rep1_time12_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.bg > JB1_WT_Rep1_time12_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.chr.bg
~/tools_av/bedGraphToBigWig JB1_WT_Rep1_time12_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_JF1_T12_JB1_WT_Rep1_str1_GRCm38.mm10.bw
~/tools_av/bedGraphToBigWig JB1_WT_Rep1_time12_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_JF1_T12_JB1_WT_Rep1_str2_GRCm38.mm10.bw
rm -r JB1_WT_Rep1_time12_AlSp_JF1_Signal.Unique.str1.out-NT.bg
rm -r JB1_WT_Rep1_time12_AlSp_JF1_Signal.Unique.str2.out-NT.bg
gzip *.bg
cd ..

cd rep2
~/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode alignReads --alignEndsType EndToEnd --runThreadN 10 --genomeDir /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-Grcm38_overlapped/ --sjdbGTFfile /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-Grcm38_overlapped/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschr.gtf --sjdbOverhang 124 --readFilesIn /home/ankitv/2019/rna-seq/time12/JB1_WT/rep2/JB1_WT_Rep2_R1_paired.fastq.gz /home/ankitv/2019/rna-seq/time12/JB1_WT/rep2/JB1_WT_Rep2_R2_paired.fastq.gz --outFileNamePrefix JB1_WT_Rep2_ --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outBAMsortingThreadN 12 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMattributes All --outSAMstrandField intronMotif --readFilesCommand gunzip -c
samtools sort -n JB1_WT_Rep2_Aligned.sortedByCoord.out.bam -o JB1_WT_Rep2_Aligned.sortedByReadname.out.bam
~/tools_av/SNPsplit_v0.3.2/SNPsplit --paired --snp_file /home/ankitv/ref_av/gencodes/gencode_M20/jf1v2_Snp-chr.GRCm38.mm10.Snpfile.txt JB1_WT_Rep2_Aligned.sortedByReadname.out.bam
samtools sort -o JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_WT_Rep2_Aligned.sortedByReadname.out.genome1.bam 
samtools sort -o JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_WT_Rep2_Aligned.sortedByReadname.out.genome2.bam
/home/ankitv/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --runThreadN 12 --inputBAMfile JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam --outWigType bedGraph --outWigNorm RPM --outWigStrand Stranded --outFileNamePrefix JB1_WT_Rep2_time12_AlSp_B6_
grep "NT" JB1_WT_Rep2_time12_AlSp_B6_Signal.Unique.str1.out.bg -v > JB1_WT_Rep2_time12_AlSp_B6_Signal.Unique.str1.out-NT.bg
sort -k1,1 -k2,2n JB1_WT_Rep2_time12_AlSp_B6_Signal.Unique.str1.out-NT.bg > JB1_WT_Rep2_time12_AlSp_B6_Signal.Unique.str1.out.-NT.sort.bg
grep "NT" JB1_WT_Rep2_time12_AlSp_B6_Signal.Unique.str2.out.bg -v > JB1_WT_Rep2_time12_AlSp_B6_Signal.Unique.str2.out-NT.bg
sort -k1,1 -k2,2n JB1_WT_Rep2_time12_AlSp_B6_Signal.Unique.str2.out-NT.bg > JB1_WT_Rep2_time12_AlSp_B6_Signal.Unique.str2.out.-NT.sort.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' JB1_WT_Rep2_time12_AlSp_B6_Signal.Unique.str1.out.-NT.sort.bg > JB1_WT_Rep2_time12_AlSp_B6_Signal.Unique.str1.out.-NT.sort.chr.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' JB1_WT_Rep2_time12_AlSp_B6_Signal.Unique.str2.out.-NT.sort.bg > JB1_WT_Rep2_time12_AlSp_B6_Signal.Unique.str2.out.-NT.sort.chr.bg
~/tools_av/bedGraphToBigWig JB1_WT_Rep2_time12_AlSp_B6_Signal.Unique.str1.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_B6_T12_JB1_WT_Rep2_str1_GRCm38.mm10.bw
~/tools_av/bedGraphToBigWig JB1_WT_Rep2_time12_AlSp_B6_Signal.Unique.str2.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_B6_T12_JB1_WT_Rep2_str2_GRCm38.mm10.bw
rm -r JB1_WT_Rep2_time12_AlSp_B6_Signal.Unique.str1.out-NT.bg
rm -r JB1_WT_Rep2_time12_AlSp_B6_Signal.Unique.str2.out-NT.bg
/home/ankitv/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --runThreadN 12 --inputBAMfile JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam --outWigType bedGraph --outWigNorm RPM --outWigStrand Stranded --outFileNamePrefix JB1_WT_Rep2_time12_AlSp_JF1_
grep "NT" JB1_WT_Rep2_time12_AlSp_JF1_Signal.Unique.str1.out.bg -v > JB1_WT_Rep2_time12_AlSp_JF1_Signal.Unique.str1.out-NT.bg
sort -k1,1 -k2,2n JB1_WT_Rep2_time12_AlSp_JF1_Signal.Unique.str1.out-NT.bg > JB1_WT_Rep2_time12_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.bg
grep "NT" JB1_WT_Rep2_time12_AlSp_JF1_Signal.Unique.str2.out.bg -v > JB1_WT_Rep2_time12_AlSp_JF1_Signal.Unique.str2.out-NT.bg
sort -k1,1 -k2,2n JB1_WT_Rep2_time12_AlSp_JF1_Signal.Unique.str2.out-NT.bg > JB1_WT_Rep2_time12_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' JB1_WT_Rep2_time12_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.bg > JB1_WT_Rep2_time12_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.chr.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' JB1_WT_Rep2_time12_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.bg > JB1_WT_Rep2_time12_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.chr.bg
~/tools_av/bedGraphToBigWig JB1_WT_Rep2_time12_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_JF1_T12_JB1_WT_Rep2_str1_GRCm38.mm10.bw
~/tools_av/bedGraphToBigWig JB1_WT_Rep2_time12_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_JF1_T12_JB1_WT_Rep2_str2_GRCm38.mm10.bw
rm -r JB1_WT_Rep2_time12_AlSp_JF1_Signal.Unique.str1.out-NT.bg
rm -r JB1_WT_Rep2_time12_AlSp_JF1_Signal.Unique.str2.out-NT.bg
gzip *.bg
cd ..
cd ..
cd JB1_ZFP57_KO
cd rep1
~/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode alignReads --alignEndsType EndToEnd --runThreadN 10 --genomeDir /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-Grcm38_overlapped/ --sjdbGTFfile /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-Grcm38_overlapped/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschr.gtf --sjdbOverhang 124 --readFilesIn /home/ankitv/2019/rna-seq/time12/JB1_ZFP57_KO/rep1/JB1_ZFP57_KO_Rep1_R1_paired.fastq.gz /home/ankitv/2019/rna-seq/time12/JB1_ZFP57_KO/rep1/JB1_ZFP57_KO_Rep1_R2_paired.fastq.gz --outFileNamePrefix JB1_ZFP57_KO_Rep1_ --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outBAMsortingThreadN 12 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMattributes All --outSAMstrandField intronMotif --readFilesCommand gunzip -c
samtools sort -n JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bam -o JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.bam
~/tools_av/SNPsplit_v0.3.2/SNPsplit --paired --snp_file /home/ankitv/ref_av/gencodes/gencode_M20/jf1v2_Snp-chr.GRCm38.mm10.Snpfile.txt JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.bam
samtools sort -o JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome1.bam 
samtools sort -o JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome2.bam
/home/ankitv/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --runThreadN 12 --inputBAMfile JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bam --outWigType bedGraph --outWigNorm RPM --outWigStrand Stranded --outFileNamePrefix JB1_ZFP57_KO_Rep1_time12_AlSp_B6_
grep "NT" JB1_ZFP57_KO_Rep1_time12_AlSp_B6_Signal.Unique.str1.out.bg -v > JB1_ZFP57_KO_Rep1_time12_AlSp_B6_Signal.Unique.str1.out-NT.bg
sort -k1,1 -k2,2n JB1_ZFP57_KO_Rep1_time12_AlSp_B6_Signal.Unique.str1.out-NT.bg > JB1_ZFP57_KO_Rep1_time12_AlSp_B6_Signal.Unique.str1.out.-NT.sort.bg
grep "NT" JB1_ZFP57_KO_Rep1_time12_AlSp_B6_Signal.Unique.str2.out.bg -v > JB1_ZFP57_KO_Rep1_time12_AlSp_B6_Signal.Unique.str2.out-NT.bg
sort -k1,1 -k2,2n JB1_ZFP57_KO_Rep1_time12_AlSp_B6_Signal.Unique.str2.out-NT.bg > JB1_ZFP57_KO_Rep1_time12_AlSp_B6_Signal.Unique.str2.out.-NT.sort.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' JB1_ZFP57_KO_Rep1_time12_AlSp_B6_Signal.Unique.str1.out.-NT.sort.bg > JB1_ZFP57_KO_Rep1_time12_AlSp_B6_Signal.Unique.str1.out.-NT.sort.chr.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' JB1_ZFP57_KO_Rep1_time12_AlSp_B6_Signal.Unique.str2.out.-NT.sort.bg > JB1_ZFP57_KO_Rep1_time12_AlSp_B6_Signal.Unique.str2.out.-NT.sort.chr.bg
~/tools_av/bedGraphToBigWig JB1_ZFP57_KO_Rep1_time12_AlSp_B6_Signal.Unique.str1.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_B6_T12_JB1_ZFP57_KO_Rep1_str1_GRCm38.mm10.bw
~/tools_av/bedGraphToBigWig JB1_ZFP57_KO_Rep1_time12_AlSp_B6_Signal.Unique.str2.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_B6_T12_JB1_ZFP57_KO_Rep1_str2_GRCm38.mm10.bw
rm -r JB1_ZFP57_KO_Rep1_time12_AlSp_B6_Signal.Unique.str1.out-NT.bg
rm -r JB1_ZFP57_KO_Rep1_time12_AlSp_B6_Signal.Unique.str2.out-NT.bg
/home/ankitv/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --runThreadN 12 --inputBAMfile JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bam --outWigType bedGraph --outWigNorm RPM --outWigStrand Stranded --outFileNamePrefix JB1_ZFP57_KO_Rep1_time12_AlSp_JF1_
grep "NT" JB1_ZFP57_KO_Rep1_time12_AlSp_JF1_Signal.Unique.str1.out.bg -v > JB1_ZFP57_KO_Rep1_time12_AlSp_JF1_Signal.Unique.str1.out-NT.bg
sort -k1,1 -k2,2n JB1_ZFP57_KO_Rep1_time12_AlSp_JF1_Signal.Unique.str1.out-NT.bg > JB1_ZFP57_KO_Rep1_time12_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.bg
grep "NT" JB1_ZFP57_KO_Rep1_time12_AlSp_JF1_Signal.Unique.str2.out.bg -v > JB1_ZFP57_KO_Rep1_time12_AlSp_JF1_Signal.Unique.str2.out-NT.bg
sort -k1,1 -k2,2n JB1_ZFP57_KO_Rep1_time12_AlSp_JF1_Signal.Unique.str2.out-NT.bg > JB1_ZFP57_KO_Rep1_time12_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' JB1_ZFP57_KO_Rep1_time12_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.bg > JB1_ZFP57_KO_Rep1_time12_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.chr.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' JB1_ZFP57_KO_Rep1_time12_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.bg > JB1_ZFP57_KO_Rep1_time12_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.chr.bg
~/tools_av/bedGraphToBigWig JB1_ZFP57_KO_Rep1_time12_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_JF1_T12_JB1_ZFP57_KO_Rep1_str1_GRCm38.mm10.bw
~/tools_av/bedGraphToBigWig JB1_ZFP57_KO_Rep1_time12_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_JF1_T12_JB1_ZFP57_KO_Rep1_str2_GRCm38.mm10.bw
rm -r JB1_ZFP57_KO_Rep1_time12_AlSp_JF1_Signal.Unique.str1.out-NT.bg
rm -r JB1_ZFP57_KO_Rep1_time12_AlSp_JF1_Signal.Unique.str2.out-NT.bg
gzip *.bg
cd ..
cd rep2
~/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode alignReads --alignEndsType EndToEnd --runThreadN 10 --genomeDir /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-Grcm38_overlapped/ --sjdbGTFfile /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-Grcm38_overlapped/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschr.gtf --sjdbOverhang 124 --readFilesIn /home/ankitv/2019/rna-seq/time12/JB1_ZFP57_KO/rep2/JB1_ZFP57_KO_Rep2_R1_paired.fastq.gz /home/ankitv/2019/rna-seq/time12/JB1_ZFP57_KO/rep2/JB1_ZFP57_KO_Rep2_R2_paired.fastq.gz --outFileNamePrefix JB1_ZFP57_KO_Rep2_ --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outBAMsortingThreadN 12 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMattributes All --outSAMstrandField intronMotif --readFilesCommand gunzip -c
samtools sort -n JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bam -o JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.bam
~/tools_av/SNPsplit_v0.3.2/SNPsplit --paired --snp_file /home/ankitv/ref_av/gencodes/gencode_M20/jf1v2_Snp-chr.GRCm38.mm10.Snpfile.txt JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.bam
samtools sort -o JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome1.bam 
samtools sort -o JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome2.bam
/home/ankitv/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --runThreadN 12 --inputBAMfile JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam --outWigType bedGraph --outWigNorm RPM --outWigStrand Stranded --outFileNamePrefix JB1_ZFP57_KO_Rep2_time12_AlSp_B6_
grep "NT" JB1_ZFP57_KO_Rep2_time12_AlSp_B6_Signal.Unique.str1.out.bg -v > JB1_ZFP57_KO_Rep2_time12_AlSp_B6_Signal.Unique.str1.out-NT.bg
sort -k1,1 -k2,2n JB1_ZFP57_KO_Rep2_time12_AlSp_B6_Signal.Unique.str1.out-NT.bg > JB1_ZFP57_KO_Rep2_time12_AlSp_B6_Signal.Unique.str1.out.-NT.sort.bg
grep "NT" JB1_ZFP57_KO_Rep2_time12_AlSp_B6_Signal.Unique.str2.out.bg -v > JB1_ZFP57_KO_Rep2_time12_AlSp_B6_Signal.Unique.str2.out-NT.bg
sort -k1,1 -k2,2n JB1_ZFP57_KO_Rep2_time12_AlSp_B6_Signal.Unique.str2.out-NT.bg > JB1_ZFP57_KO_Rep2_time12_AlSp_B6_Signal.Unique.str2.out.-NT.sort.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' JB1_ZFP57_KO_Rep2_time12_AlSp_B6_Signal.Unique.str1.out.-NT.sort.bg > JB1_ZFP57_KO_Rep2_time12_AlSp_B6_Signal.Unique.str1.out.-NT.sort.chr.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' JB1_ZFP57_KO_Rep2_time12_AlSp_B6_Signal.Unique.str2.out.-NT.sort.bg > JB1_ZFP57_KO_Rep2_time12_AlSp_B6_Signal.Unique.str2.out.-NT.sort.chr.bg
~/tools_av/bedGraphToBigWig JB1_ZFP57_KO_Rep2_time12_AlSp_B6_Signal.Unique.str1.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_B6_T12_JB1_ZFP57_KO_Rep2_str1_GRCm38.mm10.bw
~/tools_av/bedGraphToBigWig JB1_ZFP57_KO_Rep2_time12_AlSp_B6_Signal.Unique.str2.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_B6_T12_JB1_ZFP57_KO_Rep2_str2_GRCm38.mm10.bw
rm -r JB1_ZFP57_KO_Rep2_time12_AlSp_B6_Signal.Unique.str1.out-NT.bg
rm -r JB1_ZFP57_KO_Rep2_time12_AlSp_B6_Signal.Unique.str2.out-NT.bg
/home/ankitv/tools_av/STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --runThreadN 12 --inputBAMfile JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam --outWigType bedGraph --outWigNorm RPM --outWigStrand Stranded --outFileNamePrefix JB1_ZFP57_KO_Rep2_time12_AlSp_JF1_
grep "NT" JB1_ZFP57_KO_Rep2_time12_AlSp_JF1_Signal.Unique.str1.out.bg -v > JB1_ZFP57_KO_Rep2_time12_AlSp_JF1_Signal.Unique.str1.out-NT.bg
sort -k1,1 -k2,2n JB1_ZFP57_KO_Rep2_time12_AlSp_JF1_Signal.Unique.str1.out-NT.bg > JB1_ZFP57_KO_Rep2_time12_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.bg
grep "NT" JB1_ZFP57_KO_Rep2_time12_AlSp_JF1_Signal.Unique.str2.out.bg -v > JB1_ZFP57_KO_Rep2_time12_AlSp_JF1_Signal.Unique.str2.out-NT.bg
sort -k1,1 -k2,2n JB1_ZFP57_KO_Rep2_time12_AlSp_JF1_Signal.Unique.str2.out-NT.bg > JB1_ZFP57_KO_Rep2_time12_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' JB1_ZFP57_KO_Rep2_time12_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.bg > JB1_ZFP57_KO_Rep2_time12_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.chr.bg
awk '{print "chr"$1"\t"$2"\t"$3"\t"$4}' JB1_ZFP57_KO_Rep2_time12_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.bg > JB1_ZFP57_KO_Rep2_time12_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.chr.bg
~/tools_av/bedGraphToBigWig JB1_ZFP57_KO_Rep2_time12_AlSp_JF1_Signal.Unique.str1.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_JF1_T12_JB1_ZFP57_KO_Rep2_str1_GRCm38.mm10.bw
~/tools_av/bedGraphToBigWig JB1_ZFP57_KO_Rep2_time12_AlSp_JF1_Signal.Unique.str2.out.-NT.sort.chr.bg /home/ankitv/ref_av/mm10/mm10.chrom.sizes AlSp_JF1_T12_JB1_ZFP57_KO_Rep2_str2_GRCm38.mm10.bw
rm -r JB1_ZFP57_KO_Rep2_time12_AlSp_JF1_Signal.Unique.str1.out-NT.bg
rm -r JB1_ZFP57_KO_Rep2_time12_AlSp_JF1_Signal.Unique.str2.out-NT.bg
gzip *.bg
cd ..

#featurecount AS
#Note: The bam files are sorted by coordinates after SNPsplit and given as input to featurecounts. Featurecounts detect and sort/reorder automatically by readname and procees these paired-end BAM files. For details see Bulk_mm10.sh

#/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -p -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-Grcm38_overlapped/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschr.gtf -T 12 -o Alsp_JB1_time12_star-featureCounts_GRCm38.mm10.txt JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam

#Featurecount beforesplit
#/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -p -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-Grcm38_overlapped/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschr.gtf -T 12 -o Alsp_JB1_time12_star-featureCounts_GRCm38.mm10.BEFORESPLIT.txt JB1_WT_Rep1_Aligned.sortedByCoord.out.bam JB1_WT_Rep2_Aligned.sortedByCoord.out.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bam


#Igf2 matrix removed ENSMUSG00000070140.1  ENSMUSG00000115302.1 from  gtf. These genes shared a large coordinate range with Igf2 and therefore impair read assignment to Igf2, subtracted gtf file is gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschrlike-.gtf
#fgrep -f igf2likegenes.txt ./../../N-masked-JF1-GRCm38-M20-overlapped/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschr.gtf -v > gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschrlike-.gtf

#/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -p -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-GRCm38-M20-overlapped/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschrlike-.gtf -T 12 -o ./Igf2/Alsp_JB1_time12_star-featureCounts_GRCm38.mm10_Igf2.txt JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_WT_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_WT_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.B6.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByReadname.out.genome.sort.JF1.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByReadname.out.genome.sort.JF1.bam


#Featurecount beforesplit Igf2
#/home/ankitv/tools_av/subread-1.6.0-Linux-x86_64/bin/featureCounts -t exon -g gene_id -p -s 2 -a /home/ankitv/ref_av/gencodes/gencode_M20/N-masked-JF1-GRCm38-M20-overlapped/Igf2/gencode.vM20.chr_patch_hapl_scaff.annotation.chr.minuschrlike-.gtf -T 12 -o Alsp_JB1_time12_star-featureCounts_GRCm38.mm10.BEFORESPLIT_Igf2.txt JB1_WT_Rep1_Aligned.sortedByCoord.out.bam JB1_WT_Rep2_Aligned.sortedByCoord.out.bam JB1_ZFP57_KO_Rep1_Aligned.sortedByCoord.out.bam JB1_ZFP57_KO_Rep2_Aligned.sortedByCoord.out.bam

