#!/usr/bin/env python3
# Betta Stronga Fasta
# lmd: 2023-10-28

import sys
import re
import statistics
import argparse
import matplotlib.pyplot as plt
import os
import logomaker
import numpy
from Bio import AlignIO

#Functions
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
      if codon in table:
        protein+=table[codon]
      else:
        protein+='X'
  return protein

def six_frames(fasta_dict):
  frame_dict={}
  frames = [0,1,2]
  for seq_name in fasta_dict:
    frame_dict[seq_name] = {}
    for frame in frames:
      frame_dict[seq_name]["frame" + str(frame)] = fasta_dict[seq_name][frame:]
      frame_dict[seq_name]["frame_rev"+str(frame)]=reverse_complement(fasta_dict[seq_name])[frame:]
  return frame_dict

def sixframes_codon_converter(frame_dict):
  sixframes_codon_dict={}
  sixframes_codon=''
  for seq_name in frame_dict:
    sixframes_codon_dict[seq_name]={}
    for frame in frame_dict[seq_name]:
      sixframes_codon_list=[]
      for i in range(0,len(frame_dict[seq_name][frame]),3):
        if len(frame_dict[seq_name][frame][i:i+3])==3:
          sixframes_codon=frame_dict[seq_name][frame][i:i+3]
          sixframes_codon_list.append(sixframes_codon)
      sixframes_codon_dict[seq_name][frame]=sixframes_codon_list
  return (sixframes_codon_dict)

def protein_dict(sixframes_codon_dict):
  sixframes_prot_dict={}
  for seq_name in sixframes_codon_dict:
   sixframes_prot_dict[seq_name]={}
   for frame in sixframes_codon_dict[seq_name]:
     codon_str=''.join(sixframes_codon_dict[seq_name][frame])
    # print(codon_str)
     sixframes_prot_dict[seq_name][frame]=translation(codon_str)
  return sixframes_prot_dict

def get_all_CDs(sixframes_prot_dict):
  sixframes_CDs_dict={}
  for seq_name in sixframes_prot_dict:
   # sixframes_CDs_dict[seq_name]={}
    CDs_list=[]
    for frame in sixframes_prot_dict[seq_name]:
      CDs= re.findall(r"M.*?_",sixframes_prot_dict[seq_name][frame])
     # print(CDs)
      CDs_list=CDs_list+CDs
    sixframes_CDs_dict[seq_name]=CDs_list
  return sixframes_CDs_dict

def max_len_CDs(sixframes_CDs_dict):
  max_CDs_dict={}
  for seq_name in sixframes_CDs_dict:
    max_CDs=max(sixframes_CDs_dict[seq_name], key=len)
    max_CDs_dict[seq_name]=max_CDs
  return max_CDs_dict


def logo(alignment):
    alignment= AlignIO.read(open(alignment),"clustal")
    print("Alignment length %i" %alignment.get_alignment_length())
    
    seqs=list()
    
    for record in alignment:
      seqs.append(str(record.seq))
    print(str(record.seq))    
    counts_matrix=logomaker.alignment_to_matrix(seqs)
    output = logomaker.Logo(counts_matrix)
    
    #transfroming count matrix to probability matrix
    prob_mat = logomaker.transform_matrix(counts_matrix, from_type='counts' , to_type='probability')
    
    #using fade probabilityas keyword argument
    logo = logomaker.Logo(prob_mat, fade_probabilities=True, stack_order='small_on_top', font_name='Arial Rounded MT Bold')
    plt.savefig("weblogo.png")


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
    gc_content_total = round((c_count_total + g_count_total) / (c_count_total + g_count_total + a_count_total + t_count_total) *100,2) 
    print(gc_content_total)

    return gc_content_total

def GC_histo(sequences):
    # create histogram plot for GC content
    GChistogram=list()
    for sequence in sequences:
        seq = sequences[sequence].upper()
        g_count = seq.count('G')
        c_count = seq.count('C')
        a_count = seq.count('A')
        t_count = seq.count('A')
        gc_content =(c_count + g_count) / (c_count + g_count + a_count + t_count) * 100
        GChistogram.append(gc_content)
    plt.hist(GChistogram, bins=100)
    plt.savefig("gc.png")
    plt.savefig("gc.tiff")

def split(fasta_file):
    # split the file into n files (n=number of sequences)
    directory = "split_fasta_files"
    parent_dir = "."
    path = os.path.join(parent_dir, directory)
    try:
        os.mkdir(path)
        print("Directory '%s' created" %directory)
    except OSError as error:
        print("Caution to overwrite directory contents:", error.strerror)
    else:
        with open(fasta_file, "r") as file:
            i=0 
            for line in file:
                if line.startswith(">"):
                    i+=1
                    seq_name = line[1:].rstrip()    
                    outfile = open(path + "/" + "seq_{0:04d}.fasta".format(i), 'w') 
                outfile.write(line)
            outfile.close()
        print("Hurray! All files are done")

def seq_ids(fasta_dict):
    # get a list of id names
    outfile = open("seq_ids.txt", 'w')
    for seq_name in fasta_dict:
        outfile.write(f"{seq_name}\n")
    outfile.close()
    print("Seq_ids saved in 'seq_ids.txt'") 

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
    plt.xlabel("Length")
    plt.ylabel("IDs")
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
            line = line.rstrip()
            if line.startswith(">"):
                other_info = re.split(r" ", line)
                seq_name = other_info[0][1:]
            else:
                if seq_name in sequences:
                    sequences[seq_name].append(line)
                else:
                    sequences[seq_name] = list()

        for seq_name in sequences:
             sequences[seq_name] = "".join(sequences[seq_name])    
    return sequences

def split_dict(sequences):
    # split the file into n files (n=number of sequences)
    directory = "split_fasta_files"
    parent_dir = "."
    path = os.path.join(parent_dir, directory)
    try: 
        os.mkdir(path)
        print("Directory '%s' created" %directory)
    except OSError as error:
        print("Caution to overwrite directory contents:", error.strerror)
    else:
        i=0
        for seq_name in sequences:
            i+=1
            outfile = open(path + "/" + "seq_{0:04d}.fasta".format(i), "w")
            outfile.write(">" + seq_name + "\n")
            outfile.write(sequences[seq_name] + "\n")
        outfile.close()
    print("Hurray! All files are done")

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


def composition(sequences):
    #calculate aa or nt composition
    seq_keys = sequences.keys()   # gets a list of keys with seq names
    print("This is the content (%) for each sequence:")
    for seq_name in sequences:
        count_comp = {}
        for element in sequences[seq_name]:
            if element not in count_comp:
                count_comp[element] = 1 
            else:
                count_comp[element] += 1
        sorted_count = sorted(count_comp)  ### gets a list of the keys (ATGC or AAs)
    
    # prints seq name and sorted content of nt or aa    
        print(seq_name)    
        contents = []   ### gets a list of the content values
        for element in sorted_count:    
            content = round((count_comp[element]) / len(sequences[seq_name])*100, 1)    
            contents.append(content)
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
                                                                                      
######  ######  #######    #          ####### #     #    #       #   
#     # #     # #          #    #     #       #     #   # #     # #  
#     # #     # #          #    #     #       #     #  #   #   #   # 
######  ######  #####      #    #     #####   #     # #     # #     #
#       #     # #          #######    #        #   #  ####### #######
#       #     # #               #     #         # #   #     # #     #
#       ######  #               #     #######    #    #     # #     # 
"""

#MAIN

def main():
    parser = argparse.ArgumentParser(description="Thanks for using our script to parse your FASTA!", epilog=ascii_art, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", choices=["print_fasta", "all_CDs", "longest_peptide", "logo", "seq_ids", "split_dict", "split", "compositio", "art", "six_frames", "length_bar", "stats", "revcomp", "regex", "GC_content","GC_histo"], help="<-- choose one of these options to carry out on your FASTA file")
    parser.add_argument("fasta_file", help="<-- then input your FASTA")
    parser.add_argument("--pattern", help="<-- choose one of these options to carry out on your FASTA file")
    #ADD MORE ACTIONS HERE JUST LIKE ABOVE
    
    argument_input = parser.parse_args()
    

    if argument_input.action == "stats":
        fasta_input = parse_fasta(argument_input.fasta_file)
        calculate_stats(fasta_input)
    elif argument_input.action == "revcomp":
        fasta_input = parse_fasta(argument_input.fasta_file)
        print_fasta(rc_dict(fasta_input))
    elif argument_input.action=="regex":
        fasta_input = parse_fasta(argument_input.fasta_file)
        regex_seq(fasta_input,argument_input.pattern)
    elif argument_input.action=="length_bar":
        fasta_input = parse_fasta(argument_input.fasta_file)
        len_bars(fasta_input)
    elif argument_input.action=="six_frames":
        fasta_input = parse_fasta(argument_input.fasta_file)
        out = six_frames(fasta_input)
        print(out)
    elif argument_input.action=="art":
        print(ascii_art)
    elif argument_input.action=="print_fasta":
        fasta_input = parse_fasta(argument_input.fasta_file)
        print_fasta(fasta_input) 
    elif argument_input.action=="split_dict":
        fasta_input = parse_fasta(argument_input.fasta_file)
        split_dict(fasta_input)
    elif argument_input.action=="composition":
        fasta_input = parse_fasta(argument_input.fasta_file)
        composition(fasta_input)
    elif argument_input.action=="split":
        split(argument_input.fasta_file)
    elif argument_input.action=='seq_ids':
        fasta_input = parse_fasta(argument_input.fasta_file)
        seq_ids(fasta_input)
    elif argument_input.action=="GC_content":
        fasta_input = parse_fasta(argument_input.fasta_file)
        GC_content(fasta_input) 
    elif argument_input.action=="GC_histo":
        fasta_input = parse_fasta(argument_input.fasta_file)
        GC_histo(fasta_input)
    elif argument_input.action=="logo":
        logo(argument_input.fasta_file)
    elif argument_input.action=="longest_peptide":
        fasta_input = parse_fasta(argument_input.fasta_file)
        Result=six_frames(fasta_input)
        sixframes_codon_python07=sixframes_codon_converter(Result)
        prot_dict_P07=protein_dict(sixframes_codon_python07)
        CDs_dict=get_all_CDs(prot_dict_P07)
        max_len_dict=max_len_CDs(CDs_dict)
        print_fasta(max_len_dict)
    elif argument_input.action=="all_CDs":
        fasta_input = parse_fasta(argument_input.fasta_file)
        Result=six_frames(fasta_input)
        sixframes_codon_python07=sixframes_codon_converter(Result)
        prot_dict_P07=protein_dict(sixframes_codon_python07)
        CDs_dict=get_all_CDs(prot_dict_P07)
        print_fasta(CDs_dict)
 
    #IMPLEMENT ADDED ACTIONS HERE ACTIONS HERE WITH ELIF
    
if __name__ == "__main__":
    main()
