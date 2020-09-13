#!/bin/env python

# To run this code, run "python pamStats.py --ref [path to reference] --input [path to spacer bed file]" and optionally -p if prophages are present 
# only uses first fasta sequence in fasta file (if multiple)
# only analyzes first chromosome listed in bed file (if multiple)
## Make sure the following packages have been downloaded and can be imported (pandas and Biopython)
import pandas as pd
from Bio import SeqIO
import argparse

#adding reference fasta, spacer bed file, and prophage arguments
parser = argparse.ArgumentParser(description='input reference fasta and spacer bed file, output number of canonical (NGG) PAMs')
parser.add_argument('--ref', dest= "ref", help="path of reference fasta file")
parser.add_argument('--input', dest= "bed", help="path of input bed file")
parser.add_argument("-p", "--prophage", help="add prophages", action="store_true")


args = parser.parse_args()

print("")
# if prophages are present based on argument, asks for user input on phage info
if args.prophage:
	numPro = int(input ("How many prophages are in this reference? "))
	counter = 1
	currentPro = 0
	prophageNameList = []
	prophageStartList = []
	prophageEndList = []
	prophageCanonList = []
	prophageNonCanonList = []
	while counter <= numPro:
		prophageNameList.append(input("What is the name of prophage " + str(counter)+"? "))
		prophageStartList.append(input("What is the start nt position of " + prophageNameList[counter - 1]+" in the reference genome? "))
		prophageEndList.append(input("What is the end nt position of " + prophageNameList[counter - 1]+" in the reference genome? "))
		prophageCanonList.append(0)	
		prophageNonCanonList.append(0)
		
		counter += 1

counter = 0

#only pulls first FASTA record aka reference sequence from record file (not pGG32ultracr)
fastaRecord = next(SeqIO.parse(args.ref, "fasta"))

#creates string of fasta sequence
fasta = fastaRecord.seq

#imports in your BED file of reads
bed = pd.read_csv(args.bed, delimiter = "\t")


pam = str()

# Number of reads (rows) in your bed file
maxBed = len(bed.index)

# Reads
total_reads = 0
canonGenome = 0
noncanonGenome = 0
canonPro = 0
noncanonPro = 0

# If the read is in a prophage region
pro = False

#current row of the bed file
rw = 0

# starts at row(read)0 of bed file, ends when you're no longer reading lines from the first chromosome in the bed file 
# or at last bed file row
# pulls in the start/end/strand of each read (make sure your bed file matches the intended columns)
# if the start and end of a read is within a prophage region, pro designated as true
# sets current phage as the phage spacer is found in 
# if + strand then downstream pam must be gg or its noncanon, 
# otherwise must be cc upstream of start because on - strand, and fasta is top strand only
# if within a prophage region, added as canon or noncanon for prophage, 
# also adds read to separate list of individual prophage spacer counts, if multiple prophages are present

print("\n" + "Determining PAMs...")

#to only read spacers on first chromosome in bed file
firstGenome = bed.iloc[rw,0]

while (rw != maxBed) and (bed.iloc[rw,0] == firstGenome):

	total_reads +=1 
	chrStart = bed.iloc[rw,1]
	chrEnd = bed.iloc[rw, 2]
	strand = bed.iloc[rw, 5]
	
	if args.prophage:
		for i in range(numPro):
			if chrStart > int(prophageStartList[i]) and chrEnd < int(prophageEndList[i]):
				pro = True
				currentPro = i
				break
			else:
				pro = False
	
	if strand == '+':
		pam = fasta[chrEnd +1: chrEnd + 3]
		if pam.lower() == "gg":
			if pro:
				canonPro += 1
				prophageCanonList[currentPro] += 1
			else:
				canonGenome +=1
		else:
			if pro:
				noncanonPro += 1
				prophageNonCanonList[currentPro] += 1
			else:
				noncanonGenome += 1
	else:
		pam = fasta[chrStart - 3 : chrStart - 1]
		if pam.lower() == "cc":
			if pro:
				canonPro += 1
				prophageCanonList[currentPro] += 1

			else:
				canonGenome += 1
		else:
			if pro:
				noncanonPro += 1
				prophageNonCanonList[currentPro] += 1

			else:
				noncanonGenome +=1
				
	if (rw % 100000) == 0:
		print (str(rw))
	
	rw += 1

# printing our outputs!	
print("\n" + args.bed.rsplit("/", 1)[-1] + " PamStats")
print("\n Total Spacers: " + str(total_reads))
perCanTot = (float(canonGenome) + float(canonPro)) / float(total_reads) * 100
print('Percent NGG Spacers from Total= ' + "{:.2f}".format(perCanTot) + "%")

print("\n From Bacterial Genome:")
print(" - NGG Spacers: " + str(canonGenome))
print(" - Non-NGG Spacers: " + str(noncanonGenome))
perCanGen = float(canonGenome)/float(canonGenome + noncanonGenome) * 100
print(' - Percent NGG: ' + "{:.2f}".format(perCanGen) + '%')

if args.prophage:
	print("\n From Prophage Regions:")
	print(" - NGG Spacers: " + str(canonPro))
	print(" - Non-NGG Spacers: " + str(noncanonPro))
	perCanPro = float(canonPro)/float(canonPro + noncanonPro) * 100
	print(' - Percent NGG: ' + "{:.2f}".format(perCanPro) + "%")
	changingNumber = float(canonPro + noncanonPro) / float(total_reads) * 100 
	print(" - Percentage of spacers from prophages out of total: " + "{:.2f}".format(changingNumber) + "%")
	
	if numPro > 1:
		for i in range(numPro):
			print("\n From prophage " + str(prophageNameList[i]) + ":")
			print(" - NGG Spacers: " + str(prophageCanonList[i]))
			print(" - Non-NGG Spacers: " + str(prophageNonCanonList[i]))
			changingNumber = float(prophageCanonList[i] / (float(prophageCanonList[i]) + float(prophageNonCanonList[i]))) * 100
			print(" - Percent NGG: " + "{:.2f}".format(changingNumber) + "%")
			changingNumber = float(prophageCanonList[i] + prophageNonCanonList[i])/float(canonPro + noncanonPro) *100
			print(" - Percentage of prophage spacers specific for " + str(prophageNameList[i]) +": " "{:.2f}".format(changingNumber) + "%")
		
print("")