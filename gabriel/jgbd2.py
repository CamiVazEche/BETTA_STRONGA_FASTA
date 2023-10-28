#! usr/bin/env python3

import sys
import math
import re

###

#Read bionet file into dictionary of enzyme_id : cut_site
path_to_bionet_file = "/Users/jgabrielbarcia/problemsets/pfb_problemsets/repo1/python7_ps/bionet.txt"
enzymes = {}
count = 0

with open(path_to_bionet_file, "r") as bionet:
  for line in bionet:
    line = line.rstrip()
    if count < 10:
      count += 1     
      continue
    elif line == "":
      continue
    else:
      enzyme_id,cut_site = re.split(r"\s\s+", line)
      enzymes[enzyme_id] = cut_site 
      count += 1

#Get fasta file and create dictionary of gene_id : sequence
input_file = ''
genes = {}
other_info = []

try:
  input_file = sys.argv[1]
  print("User provided file name:" , input_file)
  if not re.match("^.*(\.fa|\.fasta)$", input_file):
    raise ValueError("This doesn't look like a fasta file")
  with open(input_file,"r") as fasta_file:
    for line in fasta_file:
      line = line.rstrip()
      if ">" in line:
        other_info = re.split(r" ", line)
        gene_id = other_info[0]
        genes[gene_id] = ''
      else:
        genes[gene_id] += line
except IndexError:
  print("Please provide a fasta file")
except IOError as ex:
  print("Error reading fasta file:" , input_file , ': ' , ex.strerror)


#Convert enzyme cut_sites to ACTG format (compatible with fasta file)
def RE_function(seq):
    base_dict = {'^':'^','A':'A', 'T':'T', 'C':'C', 'G':'G', 'N':'[ATCG]', 'R':'[GA]', 'Y':'[TC]', 'K':'[GT]', 'M':'[AC]', 'S':'[GC]', 'W':'[AT]', 'B':'[GTC]', 'D':'[GAT]', 'H':'[ACT]', 'V':'[GCA]'}
    return ''.join([base_dict[nt] for nt in seq])

for enzyme_id in enzymes:
  enzymes[enzyme_id] = RE_function(enzymes[enzyme_id])

#Find reverse complement of each gene to make sure all cuts are accounted for:
def reverse_complement(dna_sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna_sequence[::-1]])

genes_rev  = {gene_rev_id: reverse_complement(genes[gene_id]) for gene_rev_id, genes[gene_id] in genes.items()}
  
#Remove hats from cut_sites, replace with parentheses
for enzyme_id in enzymes:
  enzymes[enzyme_id] = re.sub("(.*)\^(.*)", r"(\1)(\2)", enzymes[enzyme_id])

#Make two new dictionaries (one for regular sequences and one for rev_compl sequences) and add cut sites in the form of hats 

results = {}
for gene_id in genes:
  results[gene_id] = {}
  for enzyme_id in enzymes:
    if re.match(r"\(", enzymes[enzyme_id]):  
      gene_hit = re.sub(enzymes[enzyme_id], r"\1^\2", genes[gene_id])
      if "^" in gene_hit:
        results[gene_id][enzyme_id] = gene_hit
        cut_seq = gene_hit.split("^")
        frag_len = [len(i) for i in cut_seq]
        print(f'{gene_id} is cut by {enzyme_id}, generating {len(cut_seq)} fragments.\nThe mean fragment length is {sum(frag_len)/len(frag_len):.2f}.\nThe minimum fragment size is {min(frag_len)} and the maximum is {max(frag_len)}\n')

results_rev = {}
for gene_rev_id in genes_rev:
  results_rev[gene_rev_id] = {}
  for enzyme_id in enzymes:
    if re.match(r"\(", enzymes[enzyme_id]):
      gene_rev_hit = re.sub(enzymes[enzyme_id], r"\1^\2", genes_rev[gene_rev_id])
      if "^" in gene_rev_hit:
        results_rev[gene_rev_id][enzyme_id] = gene_rev_hit
        cut_rev_seq = gene_rev_hit.split("^")
        frag_rev_len = [len(i) for i in cut_rev_seq]
        print(f'The reverse complement of {gene_rev_id} is cut by {enzyme_id}, generating {len(cut_rev_seq)} fragments.\nThe mean fragment length is {sum(frag_rev_len)/len(frag_rev_len):.2f}.\nThe minimum fragment size is {min(frag_rev_len)} and the maximum is {max(frag_rev_len)}\n')

