#!/bin/bash

##output total number of raw reads
if [ $1 == "count_raw" ]; then
	for i in RW*.fastq.gz;
	do 
	echo $i
	gzcat $i | awk '{s++}END{print s/4}' &>raw_counts.txt;
	done
fi

##quality trim reads	
if [ $1 == "trim" ]; then
	trimmomatic PE RW1_1.fastq.gz RW1_2.fastq.gz RW1_1.trim.fastq.gz RW1_1un.trim.fastq.gz RW1_2.trim.fastq.gz RW1_2un.trim.fastq.gz TRAILING:30 ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 MINLEN:35;
	trimmomatic PE RW2_1.fastq.gz RW2_2.fastq.gz RW2_1.trim.fastq.gz RW2_1un.trim.fastq.gz RW2_2.trim.fastq.gz RW2_2un.trim.fastq.gz TRAILING:30 ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 MINLEN:35;
	trimmomatic PE RW3_1.fastq.gz RW3_2.fastq.gz RW3_1.trim.fastq.gz RW3_1un.trim.fastq.gz RW3_2.trim.fastq.gz RW3_2un.trim.fastq.gz TRAILING:30 ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 MINLEN:35;
	trimmomatic PE RW4_1.fastq.gz RW4_2.fastq.gz RW4_1.trim.fastq.gz RW4_1un.trim.fastq.gz RW4_2.trim.fastq.gz RW4_2un.trim.fastq.gz TRAILING:30 ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 MINLEN:35;
	trimmomatic PE RW5_1.fastq.gz RW5_2.fastq.gz RW5_1.trim.fastq.gz RW5_1un.trim.fastq.gz RW5_2.trim.fastq.gz RW5_2un.trim.fastq.gz TRAILING:30 ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 MINLEN:35;
	trimmomatic PE RW6_1.fastq.gz RW6_2.fastq.gz RW6_1.trim.fastq.gz RW6_1un.trim.fastq.gz RW6_2.trim.fastq.gz RW6_2un.trim.fastq.gz TRAILING:30 ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 MINLEN:35;
	trimmomatic PE RW7_1.fastq.gz RW7_2.fastq.gz RW7_1.trim.fastq.gz RW7_1un.trim.fastq.gz RW7_2.trim.fastq.gz RW7_2un.trim.fastq.gz TRAILING:30 ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 MINLEN:35;
	trimmomatic PE RW8_1.fastq.gz RW8_2.fastq.gz RW8_1.trim.fastq.gz RW8_1un.trim.fastq.gz RW8_2.trim.fastq.gz RW8_2un.trim.fastq.gz TRAILING:30 ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 MINLEN:35;
fi

if [ $1 == "count_trim" ]; then
	for i in RW*trim.fastq.gz;
	do 
	echo $i
	gzcat $i | awk '{s++}END{print s/4}';
	done >trim_counts.txt
fi

##paired end alignments to whole genome
if [ $1 == "align" ]; then
    ~/tools/bwa/bwa mem ../ref/pGG32_NCTC8325.fa RW1_1.trim.fastq.gz RW1_2.trim.fastq.gz >../aligned/RW1-pe.sam
    ~/tools/bwa/bwa mem ../ref/pGG32_NCTC8325.fa RW2_1.trim.fastq.gz RW2_2.trim.fastq.gz >../aligned/RW2-pe.sam
    ~/tools/bwa/bwa mem ../ref/pGG32_NCTC8325.fa RW3_1.trim.fastq.gz RW3_2.trim.fastq.gz >../aligned/RW3-pe.sam
    ~/tools/bwa/bwa mem ../ref/pGG32_NCTC8325.fa RW4_1.trim.fastq.gz RW4_2.trim.fastq.gz >../aligned/RW4-pe.sam
    ~/tools/bwa/bwa mem ../ref/pGG32_NCTC8325.fa RW5_1.trim.fastq.gz RW5_2.trim.fastq.gz >../aligned/RW5-pe.sam
    ~/tools/bwa/bwa mem ../ref/pGG32_NCTC8325.fa RW6_1.trim.fastq.gz RW6_2.trim.fastq.gz >../aligned/RW6-pe.sam
    ~/tools/bwa/bwa mem ../ref/pGG32_NCTC8325.fa RW7_1.trim.fastq.gz RW7_2.trim.fastq.gz >../aligned/RW7-pe.sam
    ~/tools/bwa/bwa mem ../ref/pGG32_NCTC8325.fa RW8_1.trim.fastq.gz RW8_2.trim.fastq.gz >../aligned/RW8-pe.sam

fi

##samtools to convert alignment files to binary and index against reference
if [ $1 == "samtools" ]; then
	for i in RW*-pe.sam; 
	do 
	echo $i
	base=${i%%-*}
	echo $base
	samtools view -f 2 -F 2048 -F 4 -o $base.paired.sam $i
	awk '$2 == "83" || $2 == "163"' $base.paired.sam >$base.rev.sam
	awk '$2 == "147" || $2 == "99"' $base.paired.sam >$base.fwd.sam
	cat header.sam $base.rev.sam >$base.rev2.sam
	cat header.sam $base.fwd.sam > $base.fwd2.sam
	samtools view -b -o $base.rev.bam $base.rev2.sam
	samtools view -b -o $base.fwd.bam $base.fwd2.sam
	~/tools/samtools/samtools sort -o $base.rev.sorted.bam $base.rev.bam
	~/tools/samtools/samtools sort -o $base.fwd.sorted.bam $base.fwd.bam 
 	~/tools/samtools/samtools index $base.rev.sorted.bam
 	~/tools/samtools/samtools index $base.fwd.sorted.bam;
	done
fi

if [ $1 == "stats" ]; then
	for i in RW*sorted.bam;
	do 
	echo $i
	~/tools/samtools/samtools idxstats $i
	~/tools/samtools/samtools flagstat $i;
	done >align_counts.txt
fi

if [ $1 == "bedgraph" ]; then
	for i in *sorted.bam;
	do
	base=${i%%.*}
	j=`~/tools/samtools/samtools view $i | wc -l`
	k=`expr "1/$j*1000000" | bc -l`
	echo $j
	echo $k
	/Users/delta/tools/bedtools2/bin/bedtools genomecov -ibam $i -bg -pc -scale $k >$i.bedgraph;
	done
fi