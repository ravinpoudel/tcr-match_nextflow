#!/usr/bin/env python
# coding: utf-8



import os
import sys
import pandas as pd
from Bio import SeqIO



def get_unique_epitopes_from_tcrmatch(input_file):
    df = pd.read_csv(input_file, sep='\t')
    epitope_column = list(set(df['epitope'].tolist()))
    sublists = [element.split(',') if ',' in element else [element] for element in epitope_column]
    unique_epitopes = list(set([item for sublist in sublists for item in sublist]))
    return unique_epitopes

def generate_kmers(sequence, min_length=8, max_length=11):
    kmers = []
    for length in range(min_length, max_length + 1):
        for i in range(len(sequence) - length + 1):
            kmer = sequence[i:i + length]
            kmers.append(kmer)
    return kmers

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_tcrmatch_kmers.py <output from tcrmatch> <output_kmers_file>")
        sys.exit(1)
    input_file = sys.argv[1]
    output_file_path = sys.argv[2]
    
    unique_epitopes_from_TCRmatch = get_unique_epitopes_from_tcrmatch(input_file)

    # generate 8 - 11 bp kmers from tcrmatch results
    kmers_8_11bp = [generate_kmers(x) for x in unique_epitopes_from_TCRmatch]
    all_kmers = list(set([item for sublist in kmers_8_11bp for item in sublist]))

    with open(output_file_path, 'w') as output_file:
        for item in all_kmers:
            output_file.write(f"{item}\n")

