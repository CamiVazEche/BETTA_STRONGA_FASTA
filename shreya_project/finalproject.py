#! usr/bin/env python3

import sys
import re
import statistics
import argparse

#Functions

def parse_fasta(fasta_file):
    sequences = {}
    seq_name = ""
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith(">"):
                seq_name = line.strip()[1:]
                sequences[seq_name] = ""
            else:
                sequences[seq_name] += line.strip()
    return sequences

def calculate_stats(sequences):
    sequence_lengths = [len(seq) for seq in sequences.values()]
    total_sequences = len(sequences)
    max_length = max(sequence_lengths)
    min_length = min(sequence_lengths)
    avg_length = statistics.mean(sequence_lengths)

    sorted_lengths = sorted(sequence_lengths)
    n50 = 0
    l50 = 0
    half_total_length = sum(sequence_lengths) / 2

    for length in sorted_lengths:
        n50 += length
        if n50 >= half_total_length:
            l50 = length
            break

    print(f"Number of sequences: {total_sequences}")
    print(f"Max length: {max_length}")
    print(f"Min length: {min_length}")
    print(f"Average length: {avg_length:.2f}")
    print(f"N50: {n50}")
    print(f"L50: {l50}")
 complement_base(base):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement_dict.get(base, base)

def reverse_complement(sequences):
    for header, sequence in sequences.items():
        revcomp_seq = ''.join([complement_base(base) for base in sequence[::-1]])
        print(f"{header}\n{revcomp_seq}") 

ascii_art = """

.__  __ /\                     .__             ___.   .__  __         .__      
|__|/  |)/ ______   ___________|__| ____       \_ |__ |__|/  |_  ____ |  |__   
|  \   __\/  ___/ _/ __ \_  __ \  |/ ___\       | __ \|  \   __\/ ___\|  |  \  
|  ||  |  \___ \  \  ___/|  | \/  \  \___       | \_\ \  ||  | \  \___|   Y  \ 
|__||__| /____  >  \___  >__|  |__|\___  > /\   |___  /__||__|  \___  >___|  / 
              \/       \/              \/  )/       \/              \/     \/  




"""
table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
def translation(seq):
  protein=""
  for i in range(0,len(seq),3):
	if len(seq[i:i+3]) == 3:
		codon=seq[i:i+3]
		protein+=table[codon]
		return protein

def readingframe_1(seq):
	readingframe_1=seq[0:]
	return readingframe_1

def readingframe_2(seq):
	readingframe_2=seq[1:]
	return readingframe_2

def readingframe_3(seq):
	readingframe_3=seq[2:]
	return readingframe_3


#MAIN

def main():
    parser = argparse.ArgumentParser(description="Thanks for using our script to parse your FASTA!", epilog=ascii_art, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", choices=["stats", "revcomp"], help="<-- choose one of these options to carry out on your FASTA file")
    parser.add_argument("fasta_file", help="<-- then input your FASTA")
    #ADD MORE ACTIONS HERE JUST LIKE ABOVE
    
    argument_input = parser.parse_args()
    
    fasta_input = parse_fasta(argument_input.fasta_file)

    if argument_input.action == "stats":
        calculate_stats(fasta_input)
    elif argument_input.action == "revcomp":
        reverse_complement(fasta_input)
    #IMPLEMENT ADDED ACTIONS HERE ACTIONS HERE WITH ELIF
    
if __name__ == "__main__":
    main()
