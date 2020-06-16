'''
Created on 10 Jun 2020

@author: Clara
'''
#import the module - from '/Users/Clara/anaconda2/lib/python2.7/site-packages'
import pysam
#import the graphing software
import matplotlib
matplotlib.use('Agg') #this tells matplotlib we may not necessarily have access to a screen
from matplotlib import pyplot as plt #this is what we're going to use to plot
#define a function that we're going to call this function twice: once for plasmodium and once for cerevisiae
def getIndelFreq(fPath):
    #create two empty lists 
    insertionFreq = []
    delFreq = []    
    f = pysam.Samfile('/Users/Clara/Dropbox/barcode2little.bam','rb') #open the bam file
    for record in f:

        #uses the cigar string to calculate the number of insertions in this read
        numOfInsertions = sum(tup[1] for tup in record.cigartuples if tup[0] == 1)
        #calculates the length of the reference subsequence that this read mapped to
        lenOnRef = record.reference_length 
        #add to the insertions frequency list
        insertionFreq.append(float(numOfInsertions) / lenOnRef) 
        #same for deletions
        numOfDeletions = sum(tup[1] for tup in record.cigartuples if tup[0] == 2)  
        
        #add to deletions frequency list
        delFreq.append(float(numOfDeletions) / lenOnRef)
    f.close()
    return insertionFreq, delFreq
    print(insertionFreq)
#MAIN--------------------------
#call the function twice 
plasIns, plasDel = getIndelFreq('/Users/Clara/Dropbox/barcode2little.bam')
cerevisiaeIns, cerevisiaeDel = getIndelFreq('/Users/Clara/Dropbox/cerevisiaelittle.bam')
print(plasIns)
#do the plotting
plt.figure() #make a figure to write on
#write some plots onto the figure
plt.hist(plasIns, 50, alpha = 0.3, label = 'Falciparum Ins')
plt.hist(plasDel, 50, alpha = 0.3, label = 'Falciparum Del')
plt.hist(cerevisiaeIns, 50, alpha = 0.3, label = 'Cerevisiae Ins')
plt.hist(cerevisiaeDel, 50, alpha = 0.3, label = 'Cerevisiae Del')
#explaination of these plotting commands:
# - we're going to make a histogram out of the list in the first argument
# - the second argument specifies how many bins we want (50 in this case)
# - the third argument specifies the transparency 
# - the fourth argument specifies the label for the plot legend
#so plot housekeeping - make the legend and label the axes
plt.legend(framealpha=0.3) #make it slightly transparent so it doesn't block out any of the bars
plt.xlabel('InDels Per Reference Base')
plt.ylabel('Count')
plt.savefig('newplotOut.png') #saves the plot to a file 
plt.close()
