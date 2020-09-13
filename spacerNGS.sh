#!/bin/bash

##run scripts to demultiplex (separate out) Illumina sequencing reads based on barcodes
if [ $1 == "demux" ]; then
	python Demultiplexer.py
fi

##scripts to find spacers within sequencing reads and pull them out into a new file
if [ $1 == "parse" ]; then
	for i in 190216*; 
	do python Spacer_parser.py --input $i --bad $i.badspacers.fa --good $i.spacers.fa --meta $i.meta.txt;
	done
fi

##once we have the spacers, we want to find out where they are aligning (phage panel)	 	
if [ $1 == "align_phage" ]; then
	for i in 190216_phage_*.spacers.fa; 
	do 
	~/tools/bwa/bwa aln -n 0 ./ref/phagepanel.ref.fa $i >./align/$i.sai
	~/tools/bwa/bwa samse ./ref/phagepanel.ref.fa ./align/$i.sai $i > ./align/$i.sam
	done
fi

##alignment for prophage panel
if [ $1 == "align_newman" ]; then
	for i in 190216_prophage_1625*.spacers.fa 190216_prophage_1628*.spacers.fa; 
	do 
	~/tools/bwa/bwa aln -n 0 ./ref/prophagepanel_newman.ref.fa $i >./align/$i.sai
	~/tools/bwa/bwa samse ./ref/prophagepanel_newman.ref.fa ./align/$i.sai $i > ./align/$i.sam
	done
fi		 	

##alignment for prophage panel
if [ $1 == "align_nctc" ]; then
	for i in 190216_prophage_1626*.spacers.fa 190216_prophage_1627*.spacers.fa; 
	do 
	~/tools/bwa/bwa aln -n 0 ./ref/prophagepanel_nctc.ref.fa $i >./align/$i.sai
	~/tools/bwa/bwa samse ./ref/prophagepanel_nctc.ref.fa ./align/$i.sai $i > ./align/$i.sam
	done
fi	

##tools to view and sort alignments
if [ $1 == "samtools" ]; then
	for i in *.sam; 
	do 
	samtools view -b -o $i.bam $i
	samtools sort -o $i.sorted.bam $i.bam 
	samtools index $i.sorted.bam;
	done
fi

##output statistics for each alignment
if [ $1 == "stats" ]; then
	for i in *.sorted.bam;
	do 
	echo $i
	samtools idxstats $i
	samtools flagstat $i;
	done
fi