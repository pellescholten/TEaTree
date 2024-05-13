#!/bin/bash

# Check for the input file
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <fasta_file>"
    exit 1
fi

fasta_file=$1
output_file=$2".tsv"

# Check if the input file exists
if [ ! -f "$fasta_file" ]; then
    echo "File not found: $fasta_file"
    exit 1
fi

# Process the FASTA file
awk '/^>/ { 
        if (seqname) { 
            print seqname"\t"seqlen 
        } 
        seqname=substr($0, 2); 
        seqlen=0 
    } 
    !/^>/ { 
        seqlen+=length($0) 
    } 
    END { 
        if (seqname) print seqname"\t"seqlen 
    }' $fasta_file > $output_file
