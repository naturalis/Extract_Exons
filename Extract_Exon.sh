#!/bin/bash

# usage [genomic_fasta] [protein_ref] [output file]

genomic=$1
protein=$2
output=$3

# split the genomic file into the seperate fasta file
split -l 2 $genomic $genomic"_"

# parse through the fasta files
for fasta in $genomic"_"*; do
	VULGAR=$(exonerate --model protein2genome:bestfit $protein $fasta --exhaustive yes | grep vulgar)
	#exonerate --model protein2genome:bestfit $protein $fasta --exhaustive yes
	./parse_VULGAR.py $fasta "$VULGAR" >> $output
done

# remove redundant fasta files
rm $genomic"_"*
