#!/bin/bash
## bwa to align reads to reference genome then freebayes to call variants
## run this in your directory with the trimmed fastq files and a reference genome
echo "reference genome is the first argument"
echo "reference genome fasta should end in .fasta"
echo "the file beginning is the second argument"

bwa index $1

ls $2_*.fastq.gz | sed -e 's/\(.*\)_[R][0-9]\{1,2\}_.*$/\1_R/' | sort | uniq | while read froot ; do
    bwa mem $1 ${froot}1_paired_trimmed.fastq.gz ${froot}2_paired_trimmed.fastq.gz -t 20 |
                samtools view -b |
                samtools sort --threads 20 > ${froot}_$(basename $1 .fasta).bam
done

ls *.bam | sort | uniq | while read i ; do
       bamaddrg -b ${i} > $(basename ${i} .bam).rg.bam
done

ls *rg.bam > rgbamlist.txt

ls *.rg.bam | sort | uniq | while read i ; do
       samtools index -b -@ 20 ${i}
done

freebayes -f $1 --min-alternate-fraction 0.3 --min-coverage 10 -p 1 --bam-list rgbamlist.txt > $(basename $1 .fasta)_maf03_mincov10.vcf

##########################################################################################
#Use the following line if you are working with a reference genome that has more than one contig in it
#echo {1..32} | tr " " "\n" | parallel -j 14 "freebayes -f $1 -r contig_{} --min-alternate-fraction 0.5 --min-coverage 10 -p 1 --bam-list rgbamlist.txt > $(basename $1 .fasta)_maf05_mincov10.cont_{}.vcf"

#merge the snps from all vcf files
#ls *.vcf | uniq | while read i; do grep -v '#' $i >> $(basename $1 .fasta)_maf05_mincov10.all.vcf ; done

#check that the sample lists are in the same order for all vcfs
#cat *cont*.vcf | grep '^#C' | uniq | wc -l

#move the header into a single file, then append it to the beginning of the 'all' vcf
#grep '#' $(basename $1 .fasta)_maf05_mincov10.cont_10.vcf > header.vcf
#cat header.vcf $(basename $1 .fasta)_maf05_mincov10.all.vcf > $(basename $1 .fasta)_maf05_mincov10.all.header.vcf
