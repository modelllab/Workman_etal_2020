from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
import argparse

parser = argparse.ArgumentParser(description='input spacer aquisition sequencing reads, output parsed spacers')
parser.add_argument('--input', dest= "input", help="name of input file")
parser.add_argument('--bad', dest = "badspacers", help="what you want to name the badspacers file")
parser.add_argument('--good', dest = "output", help="what you want to name the legit spacers file")
parser.add_argument('--meta', dest = "meta", help="what you want to name the metadata output")

args = parser.parse_args()

input = open(args.input, "r")
output = open(args.output, "w")
badspacers = open(args.badspacers, "w")
meta = open(args.meta, "w")

repB = "CTCTAAAA"
repE = "GTTTTGGG"

totalreads = 0
totalspacers = 0
totalbadspacers = 0

def my_rev_complement(seq):
	return Seq(seq).reverse_complement()

for seq in input:
        totalreads += 1
        if((totalreads % 100000) == 0):
                print(totalreads)
        
        posBf = nt_search(seq, repB)
        if len(posBf) > 1:
                posEf = nt_search(seq, repE)
                if len(posEf) > 1:
                        spacer = seq[posBf[1]+9:posEf[1]]
                        spacer_rev = my_rev_complement(spacer)
                        totalspacers += 1
                        output.write(">" + str(totalspacers) + "\n" + str(spacer_rev) + "\n")
                else:
                       totalbadspacers += 1
                       badspacers.write(seq + "\n")
        else:
                totalbadspacers += 1
                badspacers.write(seq + "\n")

meta.write(args.input + "\n")                        
meta.write("total reads: " + str(totalreads) + "\n")
meta.write("total spacers: " + str(totalspacers) + "\n")
meta.write("total bad spacers: " + str(totalbadspacers) + "\n")

meta.close()
input.close()
output.close()
badspacers.close()