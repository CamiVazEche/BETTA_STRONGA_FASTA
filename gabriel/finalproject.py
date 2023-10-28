#!/usr/bin/env python3

import sys
import re
import statistics
import argparse

#Functions

#PARSE AND PRINT
def parse_fasta(fasta_file):
    sequences = {}
    seq_name = ""
    with open(fasta_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                other_info = re.split(r" ", line)
                seq_name = other_info[0][1:]
                sequences[seq_name] = ""
            else:
                sequences[seq_name] += line.strip()
    return sequences

def print_fasta(sequence_dict):
    for seq_name, sequence in sequence_dict.items():
        print(f'>{seq_name}')
        for i in range(0, len(sequence), 60):
            line = sequence[i:i + 60]
            print(line) 

#SEQUENCE STATS
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

#REVERSE COMPLEMENT - RELIES ON PRINT ABOVE
def complement_base(base):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement_dict.get(base, base)

def reverse_complement(sequences):
    revcomps_dict = {}
    for header, sequence in sequences.items():
        revcomp_seq = ''.join([complement_base(base) for base in sequence[::-1]])
        revcomps_dict[header] = revcomp_seq
    return revcomps_dict
 
        
ascii_art = """
$$$$$$$\  $$$$$$$$\ $$$$$$$\        $$\   $$\ $$$$$$$$\ $$\    $$\ $$$$$$$$\ $$$$$$$\   
$$  __$$\ $$  _____|$$  __$$\       $$ |  $$ |$$  _____|$$ |   $$ |$$  _____|$$  __$$\  
$$ |  $$ |$$ |      $$ |  $$ |      $$ |  $$ |$$ |      $$ |   $$ |$$ |      $$ |  $$ | 
$$$$$$$  |$$$$$\    $$$$$$$\ |      $$$$$$$$ |$$$$$\    \$$\  $$  |$$$$$\    $$$$$$$  | 
$$  ____/ $$  __|   $$  __$$\       \_____$$ |$$  __|    \$$\$$  / $$  __|   $$  __$$<  
$$ |      $$ |      $$ |  $$ |            $$ |$$ |        \$$$  /  $$ |      $$ |  $$ | 
$$ |      $$ |      $$$$$$$  |            $$ |$$$$$$$$\    \$  /   $$$$$$$$\ $$ |  $$ | 
\__|      \__|      \_______/             \__|\________|    \_/    \________|\__|  \__| 
                                                                                       
"""

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
        print_fasta(reverse_complement(fasta_input))
        
    #IMPLEMENT ADDED ACTIONS HERE ACTIONS HERE WITH ELIF
    
if __name__ == "__main__":
    main()
