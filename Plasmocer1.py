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
zero = 0
one = 0
two = 0
three = 0
four = 0
five = 0
six = 0
count = 0
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
        
        if Trev < 1:
            zero += 1
        elif Trev == 1:
            one +=1
        elif Trev == 2:
            two +=1
        elif Trev == 3:
            three +=1
        elif Trev == 4:
            four +=1
        elif Trev == 5:
            five +=1
        elif Trev == 6:
            six +=1
        
        percenrev = (float(Trev) / sixMerlen)
        
        tFreq.append(percenrev)
            
    elif strand == 'fwd':
        
        Tfwd = sixMerOnRef.count('T')
        Trevy.apend(Tfwd)
        percenfwd = (float(Tfwd) / sixMerlen)      
        tFreq.append(percenfwd)
        
        if Tfwd == 0:
            zero += 1
        elif Tfwd == 1:
            one +=1
        elif Tfwd == 2:
            two +=1
        elif Tfwd == 3:
            three +=1
        elif Tfwd == 4:
            four +=1
        elif Tfwd == 5:
            five +=1
        elif Tfwd == 6:
            six +=1
#print(Trevy)
Mean_Trevy=np.mean(Trevy)
print(Mean_Trevy)

Mean_Tfreq=np.mean(tFreq)
STDV_Tfreq=np.std(tFreq)
STDV_Tfreqone=np.mean(tFreq)
#print(STDV_Tfreq)
#print(Mean_Tfreq)
#print(tFreq)
#print(zero)
#print(count)
f.close() 

z=[]
z.append(float(zero)/count)
z.append(float(one)/count)
z.append(float(two)/count)
z.append(float(three)/count)
z.append(float(four)/count)
z.append(float(five)/count)
z.append(float(six)/count)
#print(z)

#do the plotting
plt.figure() #make a figure to write on
#bars = ['Zero', 'One', 'Two', 'Three', 'Four', 'Five', 'Six']
bars = [0,1,2,3,4,5,6]
plt.bar(bars, z, color =['lime', 'forestgreen', 'turquoise', 'aqua',  'blue', 'mediumslateblue'], edgecolor='navy')

#explaination of these plotting commands:
# - make histogram out of the list in the first argument
# - the second argument specifies how many bins we want (50 in this case)
# - the third argument specifies the transparency - this is useful when you have a few plots on the same axis
# - the fourth argument specifies the label for the plot legend
# plot housekeeping - make  legend and label the axes

plt.xlabel('Number of Thymidines per Sixmer')
#plt.xlim(0,0.1)
plt.ylabel('Frequency (%)')
plt.title('Thymidine Frequency Graph - Barcode 1')
plt.savefig('plasm_t2.png', dpi = 300) #saves the plot to a file 
plt.close() 