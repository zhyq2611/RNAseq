#!/usr/bin/perl -w
use strict;
my ($fq1,$fq2,$outdir,$sample,$cpu) = @ARGV;

use File::Basename;
use FindBin '$Bin';
=head1 Uasge

        perl  RNAseq_Pipeline.pl <fq1> <fq2> <outdir> <sample> <CPU number>

=cut

die `pod2text $0` if (@ARGV != 5);
system("mkdir -m 755 -p $outdir") if (!-d "$outdir");
#my $outfq1 = basename($fq1);
#my $outfq2 = basename($fq2);

######################Mapping######################################

open  OUT1,">$outdir/$sample\_map.sh";
print OUT1 "perl /haplox/users/ZhouYQ/RNA-seq/bin/correct_reads_id.pl $fq1 $fq2 $outdir/$sample\_1.fq.gz $outdir/$sample\_2.fq.gz && echo change_fq_done || exit\n";
print OUT1 "/haplox/users/hujingjing/app/hisat-0.1.6-beta/hisat --n-ceil \"L,0,0.05\" -I 50 -X 700 -t -p $cpu --known-splicesite-infile /haplox/users/hujingjing/db/hg19/Annotation/refGene.splice_site.txt --no-unal  -x /haplox/users/hujingjing/db/hg19/Index/Hisat/hg19 --novel-splicesite-outfile $outdir/$sample.novel.junc -1 $outdir/$sample\_1.fq.gz -2 $outdir/$sample\_2.fq.gz |samtools view -@ $cpu -bS - > $outdir/$sample.bam && echo map_fq_done || exit\n";
print OUT1 "samtools sort -m 80G -@ $cpu $outdir/$sample.bam -o $outdir/$sample.sort.bam && echo sort_bam_done || exit\n";
print OUT1 "samtools view -h $outdir/$sample.sort.bam |awk -F \"\\t\" '((\$0 ~ \"NH:i:1[^0-9]\"||\$0 ~ \"NH:i:1\$\") && ((and(\$2,0x2))&&(\$6==\"100M\"||(\$6~\"N\"&&\$6!~\"D\"&&\$6!~\"I\"))) )|| NF<7' > $outdir/$sample.unique.sam && echo unique_bam_done || exit\n";
print OUT1 "samtools view -@ $cpu -bS $outdir/$sample.unique.sam > $outdir/$sample.unique.bam && echo unique_sort_done || exit\n";
close OUT1;
#######################Gene expression count#########################

open  OUT2,">$outdir/$sample\_expr.sh";
print OUT2 "samtools view $outdir/$sample.unique.bam | /usr/bin/python /usr/local/bin/htseq-count -f sam -r name -i gene_id -a 0 -t exon -m intersection-nonempty -o $outdir/$sample.refgene.sam - /haplox/users/hujingjing/db/hg19/Annotation/Homo_sapiens.GRCh37_filter.gtf > $outdir/$sample.count && echo gene_count_done || exit\n";
close OUT2;

######################Novel Transcript#################################

open  OUT3,">$outdir/$sample\_noveltran.sh";
print OUT3 "/haplox/users/hujingjing/pipline/RNA-seq/bin/stringtie-1.1.1/stringtie $outdir/$sample.unique.bam -G /haplox/users/hujingjing/db/hg19/Annotation/Homo_sapiens.GRCh37_filter.gtf -j 2 -m 200 -p $cpu -o $outdir/$sample.transcript.gtf -C $outdir/$sample.transcript.txt && echo stringtie_done || exit\n";
print OUT3 "cuffcompare -r /haplox/users/hujingjing/db/hg19/Annotation/Homo_sapiens.GRCh37_filter.gtf -o $outdir/$sample\_cuffcompare $outdir/$sample.transcript.gtf && echo cuffcompare_done || exit\n";
close OUT3;

#######################Alternative Splicing############################

open  OUT4,">$outdir/$sample\_altsplicing.sh";
print OUT4 "python /haplox/users/hujingjing/project/RNA-seq/new-p/AlternativeSplice/processGTF.SAMs.py /haplox/users/hujingjing/db/hg19/Annotation/Homo_sapiens.GRCh37_filter.gtf $outdir/$sample $outdir/$sample.unique.sam $outdir && echo altsplicing_done || exit\n";
close OUT4;

#######################Call SNP Indel##################################

open  OUT5,">$outdir/$sample\_SNPIndel.sh";
print OUT5 "/usr/bin/java -Djava.io.tmpdir=/haplox/users/ZhouYQ/tmp/$sample -Xmx20g -jar /tools/GATK/picard/MarkDuplicates.jar MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=20000 INPUT=$outdir/$sample.sort.bam OUTPUT=$outdir/$sample.dedup.bam METRICS_FILE=$outdir/$sample.dedup.metrics VALIDATION_STRINGENCY=LENIENT && echo dedup_done || exit\n";
print OUT5 "/usr/bin/java -Djava.io.tmpdir=/haplox/users/ZhouYQ/tmp/$sample -jar /tools/GATK/picard/AddOrReplaceReadGroups.jar I=$outdir/$sample.dedup.bam O=$outdir/$sample.dedup_rg.bam SO=coordinate RGID=$sample RGLB=$sample RGPL=Illumina RGPU=NA RGSM=$sample && echo Add_ReadGroups_done || exit\n";
print OUT5 "samtools index $outdir/$sample.dedup_rg.bam && echo dedup_index_done || exit\n";
print OUT5 "/usr/bin/java -Djava.io.tmpdir=/haplox/users/ZhouYQ/tmp/$sample -Xmx20g -jar /tools/GATK/GenomeAnalysisTK.jar -T SplitNCigarReads -R /haplox/users/hujingjing/db/hg19/Genome/hg19.fa -I $outdir/$sample.dedup_rg.bam -o $outdir/$sample.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS && echo SplitNCigarReads_done || exit\n";
print OUT5 "/usr/bin/java -Djava.io.tmpdir=/haplox/users/ZhouYQ/tmp/$sample -Xmx20g -jar /tools/GATK/GenomeAnalysisTK.jar -R /haplox/users/hujingjing/db/hg19/Genome/hg19.fa -T BaseRecalibrator -I $outdir/$sample.split.bam -knownSites /haplox/users/ZhouYQ/Database/WES_ref/dbsnp_138.hg19.vcf -knownSites /haplox/users/ZhouYQ/Database/WES_ref/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites /haplox/users/ZhouYQ/Database/WES_ref/1000G_phase1.indels.hg19.sites.vcf -o $outdir/$sample.recal1.grp -nct $cpu && echo BaseRecalibrator_grp1_done || exit\n";
print OUT5 "/usr/bin/java -Djava.io.tmpdir=/haplox/users/ZhouYQ/tmp/$sample -Xmx20g -jar /tools/GATK/GenomeAnalysisTK.jar -R /haplox/users/hujingjing/db/hg19/Genome/hg19.fa -T BaseRecalibrator -I $outdir/$sample.split.bam -knownSites /haplox/users/ZhouYQ/Database/WES_ref/dbsnp_138.hg19.vcf -knownSites /haplox/users/ZhouYQ/Database/WES_ref/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites /haplox/users/ZhouYQ/Database/WES_ref/1000G_phase1.indels.hg19.sites.vcf -BQSR $outdir/$sample.recal1.grp -o $outdir/$sample.recal2.grp -nct $cpu && echo BaseRecalibrator_grp2_done || exit\n";
print OUT5 "/usr/bin/java -Djava.io.tmpdir=/haplox/users/ZhouYQ/tmp/$sample -Xmx20g -jar /tools/GATK/GenomeAnalysisTK.jar -R /haplox/users/hujingjing/db/hg19/Genome/hg19.fa -T PrintReads -I $outdir/$sample.split.bam -BQSR $outdir/$sample.recal2.grp -o $outdir/$sample.recal.bam -nct $cpu && echo PrintReads_done || exit\n";
print OUT5 "/usr/bin/java -Djava.io.tmpdir=/haplox/users/ZhouYQ/tmp/$sample -Xmx20g -jar /tools/GATK/GenomeAnalysisTK.jar -R /haplox/users/hujingjing/db/hg19/Genome/hg19.fa -T HaplotypeCaller -I $outdir/$sample.recal.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -nct $cpu -o $outdir/$sample.raw.vcf && echo HaplotypeCaller_done || exit\n";
print OUT5 "/usr/bin/java -Djava.io.tmpdir=/haplox/users/ZhouYQ/tmp/$sample -Xmx20g -jar /tools/GATK/GenomeAnalysisTK.jar -R /haplox/users/hujingjing/db/hg19/Genome/hg19.fa -T VariantFiltration -V $outdir/$sample.raw.vcf -window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\" -o $outdir/$sample.filter.vcf && echo VariantFiltration_done || exit\n";
print OUT5 "/thinker/net/ctDNA/annovar/table_annovar.pl $outdir/$sample.filter.vcf /thinker/net/ctDNA/annovar/humandb/ -buildver hg19 -out $outdir/$sample.filter -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all,cosmic77,clinvar_20160302 -operation g,r,r,f,f,f,f,f,f,f,f,f -nastring . -vcfinput && echo annovar_done || exit\n";
close OUT5;

#######################Fusion detection#################################
system("mkdir -m 755 -p /haplox/users/ZhouYQ/RNA-seq") if (!-d "/haplox/users/ZhouYQ/RNA-seq");
system("mkdir -m 755 -p /haplox/users/ZhouYQ/RNA-seq/fusion") if (!-d "/haplox/users/ZhouYQ/RNA-seq/fusion");
system("mkdir -m 755 -p /haplox/users/ZhouYQ/RNA-seq/fusion/fq") if (!-d "/haplox/users/ZhouYQ/RNA-seq/fusion/fq");

open  OUT6,">$outdir/$sample\_fusion.sh";
print OUT6 "perl /haplox/users/hujingjing/pipline/RNA-seq/bin/SOAPfuse-v1.26/SOAPfuse-RUN.pl -c /haplox/users/hujingjing/pipline/RNA-seq/bin/SOAPfuse-v1.26/config/config.txt -fd /haplox/users/ZhouYQ/RNA-seq/ -l $outdir/$sample\_fusion.list -o /haplox/users/ZhouYQ/RNA-seq/fusion/$sample -tp /haplox/users/ZhouYQ/RNA-seq/fusion/$sample/tmp -fm && echo SOAPfuse_done || exit\n";
close OUT6;

open  OUT6,">$outdir/$sample\_fusion.list";
print OUT6 "fusion\tfq\t$sample\t150\n";
print OUT6 "fusion\tfq\t$sample\t150\n";
close OUT6;
########################################################################

open OUT,">$outdir/$sample\_run.sh";
print OUT "#nohup sh $outdir/$sample\_fusion.sh >$outdir/$sample\_fusion.sh.o 2>$outdir/$sample\_fusion.sh.e &\n";
print OUT "nohup sh $outdir/$sample\_map.sh >$outdir/$sample\_map.sh.o 2>$outdir/$sample\_map.sh.e && echo map_pipeline_done || exit\n";
print OUT "wait\n";
print OUT "nohup sh $outdir/$sample\_expr.sh >$outdir/$sample\_expr.sh.o 2>$outdir/$sample\_expr.sh.e &\n";
print OUT "nohup sh $outdir/$sample\_noveltran.sh >$outdir/$sample\_noveltran.sh.o 2>$outdir/$sample\_noveltran.sh.e &\n";
print OUT "nohup sh $outdir/$sample\_altsplicing.sh >$outdir/$sample\_altsplicing.sh.o 2>$outdir/$sample\_altsplicing.sh.e &\n";
print OUT "nohup sh $outdir/$sample\_SNPIndel.sh >$outdir/$sample\_SNPIndel.sh.o 2>$outdir/$sample\_SNPIndel.sh.e &\n";
close OUT;

########################################################################
print "cd $outdir\nAnd\n";
print "nohup sh $outdir/$sample\_run.sh >$outdir/$sample\_run.sh.o 2>$outdir/$sample\_run.sh.e &\n";
