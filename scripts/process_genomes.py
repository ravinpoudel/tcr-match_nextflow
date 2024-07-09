#!/usr/bin/env python
# coding: utf-8



import sys
from Bio import SeqIO

def generate_kmers(sequence, min_length=8, max_length=11):
    kmers = []
    for length in range(min_length, max_length + 1):
        for i in range(len(sequence) - length + 1):
            kmer = sequence[i:i + length]
            kmers.append(kmer)
    return kmers

def generate_kmers_from_genome_protein(proteins_file, output_kmers_file):
    min_kmer_length = 8
    max_kmer_length = 11
    protein_sequences = []
    # Use SeqIO to read the proteins_file
    for record in SeqIO.parse(proteins_file, "fasta"):
        protein_sequences.append(str(record.seq))
    # Generate kmers for each protein sequence
    all_kmers = []
    for protein_sequence in protein_sequences:
        kmers = generate_kmers(protein_sequence, min_kmer_length, max_kmer_length)
        all_kmers.extend(kmers)
    # Save kmers to the specified output file
    unique_all_kmers = list(set(all_kmers))
    with open(output_kmers_file, 'w') as file:
        for kmer in unique_all_kmers:
            file.write(f"{kmer}\n")
    print(f"Total number of elements in all_kmers: {len(all_kmers)}")
    print(f"Total number of elements in unique_all_kmers: {len(unique_all_kmers)}")
    print(f"All unique kmers saved to {output_kmers_file}")
    return unique_all_kmers

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <proteins_file> <output_kmers_file>")
        sys.exit(1)
    proteins_file = sys.argv[1]
    output_kmers_file = sys.argv[2]
    generate_kmers_from_genome_protein(proteins_file, output_kmers_file)
