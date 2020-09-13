#!/bin/bash

##download dsk from precompiled binary here (dskv-2.3.0-bin-Darwin.tar.gz): https://github.com/GATB/dsk/releases
##to run this code, type in 'bash dsk.sh'. Update paths to file and dsk program as needed

i = "/path/to/fasta/name.fa"

./dsk -kmer-size 8 -abundance-min 2 -file $i.fa 
./dsk2ascii -file $i.h5 -out kmerlist
less kmerlist


