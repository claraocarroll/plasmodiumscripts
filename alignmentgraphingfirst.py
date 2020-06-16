'''
Created on 15 Jun 2020

@author: Clara
'''
f = open('route/to/alignment/file', 'r')

import matplotlib
matplotlib.use('Agg') # may not necessarily have access to a screen
from matplotlib import pyplot as plt # to plot

#read in header (first line)
#print(f.readline().strip())
#create empty list
prev = []
insertionFreq = []
delFreq = []
insCount = 0
delCount = 0

for line in f:
        if line[0] == '#':
            continue
        
        if line.startswith('>'):
            splitline = line.rstrip().split(' ')
            readID = splitline[0][1:]
            chromosome = splitline[1]
            refStart = int(splitline[2])
            refEnd = int(splitline[3])
            strand = splitline[4]
            numBase = refEnd - refStart
            prevPos = -1
            prev = []
            
        else:
            splitline = line.rstrip().split('\t')
            posOnRef = int(splitline[0])
            #posOnRef.append(prev)
            prev.append(posOnRef)
            sixMerOnRef = splitline[1]
            nanoMagn = float(splitline[2])
            nanoSigLen = float(splitline[3])
            sixMerOnStrand = splitline[4]
            levelMean = float(splitline[5])
            levelStDev = float(splitline[5])
            pos = int(splitline[0])
          
            #search for insertions
            if sixMerOnStrand == 'NNNNNN': 
                insCount += 1
             
            #search for deletions
            if prevPos != -1:
                gap = abs(pos - prevPos)
                if gap > 1:
                    x = gap - 1
                    delCount += x
                    delFreq.append(float(delCount) / numBase)
            
            if numBase != -1:
                insertionFreq.append(float(insCount) / numBase)
             #reset counters   
            insCount = 0
            delCount = 0 
            delFreq.append(float(delCount) / numBase)
            #reset position
            prevPos = pos
            
f.close()

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
plt.ylabel('Count')
plt.savefig('ins_3.pdf') #saves the plot to a file 
plt.close() 
