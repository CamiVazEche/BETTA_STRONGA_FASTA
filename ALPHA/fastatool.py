#!/usr/bin/env python3
# Betta Stronga Fasta
# lmd: 2023-10-28

import sys
import re
import statistics
import argparse
import matplotlib.pyplot as plt

#Functions
def GC_content(sequences):
        # calculate gc statistics for a FASTA dictionary
    g_count_total = 0 
    c_count_total = 0 
    a_count_total = 0 
    t_count_total = 0 
    for sequence in sequences:
        seq = sequences[sequence].upper()
        g_count_total += seq.count('G')
        c_count_total += seq.count('C')
        a_count_total += seq.count('A')
        t_count_total += seq.count('T')
    gc_content_total = (c_count_total + g_count_total) / (c_count_total + g_count_total + a_count_total + t_count_total)
    return gc_content_total

def regex_seq(sequences, pattern):
    # search for regex in sequence
    print("\t".join(["Name", "Start site","Stop site"]))
    for header in sequences:
        for match_obj in re.finditer(pattern, sequences[header]):  
            print("\t".join([header,str(match_obj.start()),str( match_obj.end())])) 

def len_bars(fasta_input):
    # create bar plot of sequence lengths
    lengths={}
    for header in fasta_input:
        x = header
        y= len(fasta_input[header])
        plt.barh(x,y)
    plt.xlabel("ID")
    plt.ylabel("Nucleotides")
    plt.title("Sequence lengths")
    plt.savefig("bar.png")
    
def six_frames(fasta_dict):
  frame_dict={}
  frames = [0,1,2]
  for seq_name in fasta_dict:
    frame_dict[seq_name] = {}
    for frame in frames:
      frame_dict[seq_name]["frame" + str(frame)] = fasta_dict[seq_name][frame:]
      frame_dict[seq_name]["frame_rev"+str(frame)]=reverse_complement(fasta_dict[seq_name])[frame:]
  return frame_dict

def reverse_complement(sequence):
  revcomp_seq = ''.join([complement_base(base) for base in sequence[::-1]])
  return(revcomp_seq)

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
    gc = GC_content(sequences)

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
    print(f'GC: {gc:.4%}')


def composition(sequences):
    #calculate aa or nt composition
    for seq_name in sequences:
        count_comp = {}
        perc_comp = {}
        for element in sequences[seq_name]:
            if element not in count_comp:
                count_comp[element] = 1 
            else:
                count_comp[element] += 1
    
        for element in count_comp:    
            content = (count_comp[element]) / len(sequences[seq_name])*100    
            print(element, " ", content)

#REVERSE COMPLEMENT - RELIES ON PRINT ABOVE
def complement_base(base):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement_dict.get(base, base)

def rc_dict(sequences):
    revcomps_dict = {}
    for header in sequences:
        revcomps_dict[header] = reverse_complement(sequences[header])
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
    parser.add_argument("action", choices=["composition", "art", "six_frames", "length_bar", "stats", "revcomp", "regex"], help="<-- choose one of these options to carry out on your FASTA file")
    parser.add_argument("fasta_file", help="<-- then input your FASTA")
    parser.add_argument("--pattern", help="<-- choose one of these options to carry out on your FASTA file")
    #ADD MORE ACTIONS HERE JUST LIKE ABOVE
    
    argument_input = parser.parse_args()
    
    fasta_input = parse_fasta(argument_input.fasta_file)

    if argument_input.action == "stats":
        calculate_stats(fasta_input)
    elif argument_input.action == "revcomp":
        print_fasta(rc_dict(fasta_input))
    elif argument_input.action=="regex":
        regex_seq(fasta_input,argument_input.pattern)
    elif argument_input.action=="length_bar":
        len_bars(fasta_input)
    elif argument_input.action=="six_frames":
        out = six_frames(fasta_input)
        print(out)
    elif argument_input.action=="art":
        print(ascii_art)
    elif argument_input.action=="composition":
        composition(fasta_input)
        
    #IMPLEMENT ADDED ACTIONS HERE ACTIONS HERE WITH ELIF
    
if __name__ == "__main__":
    main()
