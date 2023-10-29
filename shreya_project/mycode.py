#!/usr/bin/env python3
import re
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
#print(aa_seq)
  

def parse_fasta(fasta_file):
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
#print(output)

    



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
#print(Result)

'''def create_CDs(fasta_dict):
  CDs_dict={}
  for seq_name in frame_dict:
    for frame in frame_dict[seq_name]:
      frame_dict[seq_name][frame]=list(sequence)
      print(frame_dict)
    start_codon=re.search(r"(.ATG)",fasta_dict[seq_name])
    CDs_start=fasta_dict[seq_name].index('ATG')
    CDs=fasta_dict[seq_name][CDs_start:]
    CDs_dict[seq_name]=CDs
  return CDs_dict
CDs_rec=create_CDs(output)
print(CDs_rec)'''
def codon_convert(fasta_dict):
  codon_dict={}
  codon_list=[]
  codon=''
  for seq_name in fasta_dict:
    for i in range (0,len(fasta_dict[seq_name]),3):
      if len(fasta_dict[seq_name][i:i+3]) == 3:
        codon=fasta_dict[seq_name][i:i+3]      
        codon_list.append(codon)
        codon_dict[seq_name]=codon_list
  return(codon_dict)
codon_fasta07=codon_convert(output)
#print(codon_fasta07)

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
sixframes_codon_python07=sixframes_codon_converter(Result)
print(sixframes_codon_python07)

'''def create_CDs_dict(codon_dict):
  CDs_dict={}
  
  for seq_name in codon_dict:
    start_codon=re.search(r"(.ATG)",codon_dict[seq_name])
    CDs_start=codon_dict[seq_name].index('ATG')
    CDs=codon_dict[seq_name][CDs_start:]     
    CDs_dict[seq_name]=CDs
  return CDs_dict
CDs_rec=create_CDs_dict(codon_fasta07)
print(CDs_rec)'''

'''def create_CDs_dict (sixframes_codon_dict):
  sixframes_CDs_dict={}
  for seq_name in sixframes_codon_dict:   
    sixframes_CDs_dict[seq_name]={}
    for frame in sixframes_codon_dict[seq_name]:
     CDs_start=sixframes_codon_dict[seq_name][frame].index('ATG')
     CDs=sixframes_codon_dict[seq_name][frame][CDs_start:]
     sixframes_CDs_dict[seq_name][frame]=CDs
  return sixframes_CDs_dict
CDs_dict_pythonsixf=create_CDs_dict(sixframes_codon_python07)
print(CDs_dict_pythonsixf)'''

def protein_dict(sixframes_codon_dict):
  sixframes_prot_dict={}
  for seq_name in sixframes_codon_dict:
   sixframes_prot_dict[seq_name]={}
   for frame in sixframes_codon_dict[seq_name]:
     codon_str=''.join(sixframes_codon_dict[seq_name][frame])
    # print(codon_str)
     sixframes_prot_dict[seq_name][frame]=translation(codon_str)
  return sixframes_prot_dict
prot_dict_P07=protein_dict(sixframes_codon_python07)
print(prot_dict_P07)

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
CDs_dict=get_all_CDs(prot_dict_P07)
print(CDs_dict)
'''def rm_small_CDs(sixframes_CDs_dict):
  sixframes_CDs_dict_rms={}
  small_CDs_list=[]
  avg_CDs_list=[]
  for seq_name in sixframes_CDs_dict:
    if len(sixframes_CDs_dict[seq_name])<5:
       small_CDs=(sixframes_CDs_dict[seq_name])
       small_CDs_list.append(small_CDs)
       avg_CDs=sixframes_CDs_dict[seq_name].remove(smallest_CDs)
       avg_CDs_list.append(avg_CDs)
       sixframes_CDs_dict_rms[seq_name]=avg_CDs_list
  return sixframes_CDs_dict_rms
rm_CDs_P07=rm_small_CDs(CDs_dict)
print(rm_CDs_P07)'''
def max_len_CDs(sixframes_CDs_dict):
  max_CDs_dict={}
  for seq_name in sixframes_CDs_dict:
    max_CDs=max(sixframes_CDs_dict[seq_name], key=len)    
    max_CDs_dict[seq_name]=max_CDs
  return max_CDs_dict
max_len_dict=max_len_CDs(CDs_dict)
print(max_len_dict)
