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
# define another function 
def getIndelFreq(fPath):
    #create empty list
    plasDel = []    
    f = pysam.Samfile('/Users/Clara/Dropbox/barcode2little.bam','rb') #open the bam file
    for record in f:
        #calculates the length of the reference subsequence that this read mapped to
        lenOnRef = record.reference_length 
        #same for deletions
        numOfDeletions = sum(tup[1] for tup in record.cigartuples if tup[0] == 2)
        #add to deletions frequency list
        plasDel.append(float(numOfDeletions) / lenOnRef)
        #print(plasIns)
    f.close()
    return plasDel
    
#MAIN--------------------------
#call the function 
plasDel = getIndelFreq('/Users/Clara/Dropbox/barcode2little.bam')

#define the function 
def getIndellFreqq(fPath):
    #create empty lists 
    cerevisiaeDel = []    
    f = pysam.Samfile('/Users/Clara/Dropbox/cerevisiaelittle.bam','rb') #open the bam file
    for record in f:
        #write something that calculates the length of the reference subsequence that this read mapped to
        lenOnRef = record.reference_length 
        # for deletions
        numOfDeletions = sum(tup[1] for tup in record.cigartuples if tup[0] == 2)
        #add to deletions frequency list
        cerevisiaeDel.append(float(numOfDeletions) / lenOnRef)
    f.close()
    return cerevisiaeDel

cerevisiaeDel = getIndellFreqq('/Users/Clara/Dropbox/cerevisiaelittle.bam')
#do the plotting
plt.figure() #make a figure to write on
#write some plots onto the figure
plt.hist(plasDel, 50, alpha = 0.3, label = 'Falciparum Del')
plt.hist(cerevisiaeDel, 50, alpha = 0.3, label = 'Cerevisiae Del')
#explaination of these plotting commands:
# -  make a histogram out of the list in the first argument
# - the second argument specifies how many bins we want (50 in this case)
# - the third argument specifies the transparency -  useful when you have a few plots on the same axis
# - the fourth argument specifies the label for the plot legend
#so plot housekeeping - make the legend and label the axes
plt.legend(framealpha=0.3) #make it slightly transparent so it doesn't block out any of the bars
plt.xlabel('InDels Per Reference Base')
plt.ylabel('Count')
plt.savefig('deletionsplotOut.png') #saves the plot to a file 
plt.close() 
