#!/usr/bin/env python3
# Betta Stronga Fasta



import sys
import re
import statistics
import argparse
import matplotlib.pyplot as plt
#import biopython
#from weblogo import *
import pandas
import logomaker
import numpy
test_file=sys.argv[1]

from Bio import AlignIO

alignment= AlignIO.read(open(test_file),"clustal")
print("Alignment length %i" %alignment.get_alignment_length())

seqs=list()

for record in alignment:
 print(record.seq + "" + record.id)
 seqs.append(str(record.seq))

print(seqs)
counts_matrix=logomaker.alignment_to_matrix(seqs)
print(counts_matrix)
 #counts_mat.head()
#count matrix to probability matrix
output = logomaker.Logo(counts_matrix)



#transfroming count matrix to probability matrix
prob_mat = logomaker.transform_matrix(counts_matrix, from_type='counts' , to_type='probability')
#print(prob_mat)

#logo = logomaker.Logo(prob_mat)
#logo = logomaker.Logo(prob_mat)

#using fade probabilityas keyword argument
logo = logomaker.Logo(prob_mat, fade_probabilities=True, stack_order='small_on_top', font_name='Arial Rounded MT Bold')
plt.savefig("weblogo.png")
#plt.show()

# seqList = weblogo.seq_io.read(afile)



