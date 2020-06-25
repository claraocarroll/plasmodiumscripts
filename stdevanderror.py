'''
Created on 25 Jun 2020

@author: Clara
'''

f = open('/path/to/file', 'r')
#import the graphing software
import numpy as np

import matplotlib

matplotlib.use('Agg') # may not necessarily have access to a screen
from matplotlib import pyplot as plt # to plot
#read in header (first line)
#create empty list
tFreq = []
zero = 0
one = 0
two = 0
three = 0
four = 0
five = 0
six = 0
count = 0
countrr = 0
probtot0 = 0
probtot1 = 0
probtot2 = 0
probtot3 = 0
probtot4 = 0
probtot5 = 0
probtot6 = 0
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
    probBrdu = float(splitline[1])
    sixMerOnRef = splitline[2]
    sixMerlen = len(sixMerOnRef)
    count +=6
    
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
        motif = 'T'
        
        if motif in rev_comp:
            countrr=0
        else:
            zero += 1
            probtot0 += probBrdu
                       
        if Trev == 1:
            one +=1
            probtot1 += probBrdu
        elif Trev == 2:
            two +=1
            probtot2 += probBrdu
        elif Trev == 3:
            countrr+=1
            probtot3 += probBrdu
            three +=1
        elif Trev == 4:
            probtot4 += probBrdu
            four +=1
        elif Trev == 5:
            probtot5 += probBrdu
            five +=1
        elif Trev == 6:
            probtot6 += probBrdu
            six +=1
        
        percenrev = (float(Trev) / sixMerlen)
        
        tFreq.append(percenrev)
            
    elif strand == 'fwd':
        
        Tfwd = sixMerOnRef.count('T')
        Trevy.append(Tfwd)
        percenfwd = (float(Tfwd) / sixMerlen)      
        tFreq.append(percenfwd)     
        motif = 'T'
        
        if motif in sixMerOnRef:
            countrr = 0

        else:
            zero += 1
            probtot0 +=probBrdu   
        
        if Tfwd == 1:
            one +=1
            probtot1 +=probBrdu
        elif Tfwd == 2:
            two +=1
            probtot2 +=probBrdu
        elif Tfwd == 3:
            three +=1
            probtot3 +=probBrdu
        elif Tfwd == 4:
            probtot4 +=probBrdu
            four +=1
        elif Tfwd == 5:
            probtot5 +=probBrdu
            five +=1
        elif Tfwd == 6:
            probtot6 +=probBrdu
            six +=1
#print(Trevy)
Mean_Trevy=np.mean(Trevy)
Mean_Tfreq=np.mean(tFreq)
STDV_Tfreq=np.std(tFreq)
STDV_Tfreqone=np.mean(tFreq)
f.close() 

#counter that keeps track of how many three-thymidine sixmers we have
# and then add probBrdU to a running total of probabilities that we have for three-thymidine sixmers.
#Then at the end, take (running total of probBrdU for three-thymidine sixMers) / (total number of three-thymidine sixMers)

z=[]
#z.append(float(probtot0)/zero)
z.append(float(probtot1)/one)
z.append(float(probtot2)/two)
z.append(float(probtot3)/three)
z.append(float(probtot4)/four)
z.append(float(probtot5)/five)
z.append(float(probtot6)/six)
print(z)
STDV_z=np.std(z)
print(STDV_z)
bars = [1,2,3,4,5,6]

error = [STDV_z]

fig, ax = plt.subplots()
# to change from a line plot with error bars to a bar plot with error bars remove the word 'error' from the below statement
ax.errorbar(bars, z, yerr=STDV_z, alpha=0.5, ecolor='black' )
ax.set_ylabel('Frequency')
ax.set_xticks(bars)
ax.set_title('Thymidine Frequency Graph - Plasmodium')
ax.yaxis.grid(True)

# Save the figure and show
plt.tight_layout()
plt.savefig('Plasmerror1.png') 
plt.close() 
