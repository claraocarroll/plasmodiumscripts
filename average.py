'''
Created on 24 Jun 2020

@author: Clara
'''
f = open('/Users/Clara/Dropbox/hpcoca1/output1_1000.detect', 'r')
#f = open('/rds/project/mb915/rds-mb915-data_transfer/test.align', 'r') 
#import the graphing software
import numpy as np
import matplotlib
matplotlib.use('Agg') # may not necessarily have access to a screen
from matplotlib import pyplot as plt # to plot
#read in header (first line)
#create empty list
tFreq = []
Trevy = []
for line in f:
  if line[0] == '#':
        continue
  if line.startswith('>'):
        #wrap up from the read before, if there is one
        prevPos = -1
        splitline = line.rstrip().split(' ')
        readID = splitline[0][1:]
        chromosome = splitline[1]
        refStart = int(splitline[2])
        refEnd = int(splitline[3])
        strand = splitline[4]
        
  else:
    splitline = line.rstrip().split('\t')
    posOnRef = int(splitline[0])
    numB = float(splitline[1])
    sixMerOnRef = splitline[2]
    sixMerlen = len(sixMerOnRef)
    
    if strand == 'rev':
        # reverse complement
        # dictionary for complement nucleotides
        complement = { 'A' : 'T', 'C' : 'G', 'G' : 'C', 'T' : 'A' }

        # initialise string variable to hold reverse complement of dna string
        rev_comp = ''

        # reverse sixmer string
        SixMer_rev = sixMerOnRef[::-1]
        
        # cycle through  string and form complement
        for nuc in SixMer_rev :
            # add complement of nucleotide to sequence string
            #rev_comp += complement[nuc]
            rev_comp += complement[nuc]
           
        # get number of T's and store in variable
        Trev = rev_comp.count('T')
        Trevy.append(Trev)
        
        
    elif strand == 'fwd':
        Tfwd = sixMerOnRef.count('T')
        Trevy.append(Tfwd)
       
        
        
#print(Trevy)
Mean_Trevy=np.mean(Trevy)
print(Mean_Trevy)


#print(STDV_Tfreq)
#print(Mean_Tfreq)
#print(tFreq)
#print(zero)
#print(count)
f.close()
