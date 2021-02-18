import pandas as pd
import os
import BitVector as bv

drivers = pd.read_csv('TableS3_panorama_driver_mutations_ICGC_samples.controlled.tsv', delimiter = '\t')

f = open("K562_distal_both_FDR_0.1.txt", "r")

def essentialElementReadIn(self, input):
    #idea: read in file, split by chromosome, parse into 2D array of strings, where each subarray contains the coordinates of each chromosome
    return #2D array of strings

def constructBitVectorArray(self, input):
    #idea: create one bit vector for every chromosome, vector should be length 24, with 0 as a dummy element
    #run through call helper
    return #array of bit vector objects

def constructBitVector(self, input):
    #helper method
    #idea: create bit vector initialized w/ 0's, run through files and initialize all values within coordinate range to 1's
    #again, 0 is dummy value to minimize indexing confusion
    return #a bit vector

def writeFiles(self, input):
    #write out the intersection based on bit vector??
    #returns a text file

#bad file read code
#filename = 'K562_distal_depleted_FDR_0.1' # select file by type file name
#with open(filename) as file_obj: # open file with open() function and make alias name from filename (file_obj)
    #for line in file_obj: # using for looping to read entire content line-by-line
        #print(line.rstrip())
        
#potential write file code:
