# Demultiplex internal indexes
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search

input = open("RW1_S1_L001_R1_001.fastq","r")
output1 = open('190216_prophage_1625minus', 'w')
output2 = open('190216_prophage_1625plus', 'w')
output3 = open('190216_jw_10aminus', 'w')
output4 = open('190216_phage_116', 'w')
output5 = open('190216_phage_68', 'w')
output6 = open('190216_phage_99', 'w')


records = SeqIO.parse(input,'fastq')

bar1 = 'ATCA'
bar2 = 'TTAG'
bar3 = 'ACTT'
bar4 = 'GATC'
bar5 = 'TAGC'
bar6 = 'GGCT'

total_reads = 0
output1_reads = 0
output2_reads = 0
output3_reads = 0
output4_reads = 0
output5_reads = 0
output6_reads = 0

for r in records:
        total_reads += 1
        seq = str(r.seq)
        barcode = seq[0:4]
        if barcode == bar1:
                output1_reads += 1
                output1.write(">" + r.id + '\n' + seq + '\n')
        if barcode == bar2:
                output2_reads += 1
                output2.write(">" + r.id + '\n' + seq + '\n')
        if barcode == bar3:
                output3_reads += 1
                output3.write(">" + r.id + '\n' + seq + '\n')
        if barcode == bar4:
                output4_reads += 1
                output4.write(">" + r.id + '\n' + seq + '\n')
        if barcode == bar5:
                output5_reads += 1
                output5.write(">" + r.id + '\n' + seq + '\n')
        if barcode == bar6:
                output6_reads += 1
                output6.write(">" + r.id + '\n' + seq + '\n')        

print(total_reads)
print(output1_reads)
print(output2_reads)
print(output3_reads)
print(output4_reads)
print(output5_reads)
print(output6_reads)
                                
input.close()
output1.close()
output2.close()
output3.close()
output4.close()
output5.close()
output6.close()