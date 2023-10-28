#!/usr/bin/env python3
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
seq='ATGCGTGACGTTTGACGGTAGGCGTGAGGCGATTTAGCGGTAGGACGATGCGATGATGATG'
aa_seq= translation(seq)
  

'''def parse_fasta(fasta_file):
  sequences={}
  with open(fasta_file, "r") as file_obj:
    for line in file_obj:
      if line.startswith('>'):
        seq_name=line
        sequences[seq_name]=''# also, we can make a list instead of the string
      if not line.startswith('>'):
        line=line.rstrip()
        sequence=line
        sequences[seq_name] += sequence
                                       #we can join the list here to get a string
  return sequences
output=parse_fasta("Python_07.fasta")
print(output)'''

    



def six_frames(fasta_dict):
  frame_dict={}
  frames = [0,1,2]
  for seq_name in fasta_dict:
    frame_dict[seq_name] = {}
    for frame in frames:
      frame_dict[seq_name]["frame" + str(frame)] = fasta_dict[seq_name][frame:]
      frame_dict[seq_name]["frame_rev"+str(frame)]=reverse_complement(fasta_dict[seq_name])[frame:]
  return frame_dict

def complement_base(base):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement_dict.get(base, base)

def reverse_complement(sequence):
  revcomp_seq = ''.join([complement_base(base) for base in sequence[::-1]])
  return(revcomp_seq)

Result=six_frames(output)
print(Result)
