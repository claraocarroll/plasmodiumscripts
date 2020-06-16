'''
Created on 15 Jun 2020

@author: Clara
'''
f = open('/home/hpcoca1/rds/hpc-work/resultsalign_2020_06_08/output7.align', 'r')
#f = open('/rds/project/mb915/rds-mb915-data_transfer/test.align', 'r') 
#import the graphing software

import matplotlib
matplotlib.use('Agg') # may not necessarily have access to a screen
from matplotlib import pyplot as plt # to plot
#read in header (first line)
#create empty lists
insertionFreq = []
delFreq = []
insCount = 0
delCount = 0
numBase = -1
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
    if numBase != -1:
        delFreq.append(float(delCount) / numBase)
        insertionFreq.append(float(insCount) / numBase)
    numBase = refEnd - refStart 
    
    #reset counters for the next read
    
    delCount = 0
    insCount = 0
    
  else:
    splitline = line.rstrip().split('\t')
    posOnRef = int(splitline[0])
    sixMerOnRef = splitline[1]
    nanoMagn = float(splitline[2])
    nanoSigLen = float(splitline[3])
    sixMerOnStrand = splitline[4]
    levelMean = float(splitline[5])
    levelStDev = float(splitline[5])
    pos = int(splitline[0])
    
    if prevPos != -1:
      gap = abs(pos - prevPos)
      if gap > 1:
        delCount += gap - 1
    
    prevPos = pos
    
    #search for insertions
    if sixMerOnStrand == 'NNNNNN':
      insCount += 1
      
f.close()
print(insertionFreq[0:10], 'insertion')
print(delFreq[0:10], 'deletion')
   
print(len(insertionFreq), 'insertion length')
print(len(delFreq), 'deletion length')
#do the plotting
plt.figure() #make a figure to write on
#write some plots onto the figure
plt.hist(insertionFreq, 50, alpha = 0.3, label = 'Falciparum Ins')
plt.hist(delFreq, 50, alpha = 0.3, label = 'Falciparum Del')
#explaination of these plotting commands:
# - make histogram out of the list in the first argument
# - the second argument specifies how many bins we want (50 in this case)
# - the third argument specifies the transparency - this is useful when you have a few plots on the same axis
# - the fourth argument specifies the label for the plot legend
# plot housekeeping - make  legend and label the axes
plt.legend(framealpha=0.3) #make it slightly transparent so it doesn't block out any of the bars
plt.xlabel('InDels per Reference Base')
plt.xlim(0,0.1)
plt.ylabel('Count')
plt.savefig('outboemo9_15.pdf') #saves the plot to a file 
plt.close() 
