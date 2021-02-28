import pandas as pd
import os
from BitVector import BitVector
import sys
import time 


def essentialElementReadIn(file_path):
    df = pd.read_csv(file_path, sep = "\t", names = ['Chr_Num', 'Start', 'End'])

    #cleaning data
    df.Chr_Num = df.Chr_Num.str[3:] #getting rid of Chr tag
    df[['Start', 'End']] = df[["Start", 'End']].apply(pd.to_numeric) #changing start and end coords to ints

    #creating map
    chromosome_coord_map = {}

    for index, row in df.iterrows():
        chrom_num = row['Chr_Num']
    
        if chrom_num not in chromosome_coord_map.keys():
            chromosome_coord_map[chrom_num] = []
        chromosome_coord_map[chrom_num].extend([row['Start'], row['End']])

    return chromosome_coord_map

def chromosomeLengthReadIn(file_path):
    df = pd.read_csv(file_path, names = ['Chr_Num', 'Length'])
    df['Length'] = pd.to_numeric(df['Length']) 

    chromosome_length_map = {}

    for index, row in df.iterrows():
        chrom_num = row['Chr_Num']
        chromosome_length_map[chrom_num] = row['Length']

    return chromosome_length_map


def constructBitVector(chr_num):
#     #helper method
#     #idea: create bit vector initialized w/ 0's, run through files and initialize all values within coordinate range to 1's
#     #again, 0 is dummy value to minimize indexing confusion
#     return a bit vector
    chromosome_coord_map = essentialElementReadIn("data/K562_distal_both_FDR_0.1.txt")

    chromosome_length_map = chromosomeLengthReadIn("data/ChromosomeLengths.txt")

    currIndex = 0

    bv  = BitVector(size = chromosome_length_map[chr_num] + 1)

    coord_lst = chromosome_coord_map[chr_num]

    for i in range(chromosome_length_map[chr_num]):
        if (i < coord_lst[currIndex] and currIndex % 2 == 0):
            bv[i + 1] = 0
        elif(i < coord_lst[currIndex] and currIndex % 2 == 1):
            bv[i + 1] = 1
        elif(currIndex + 1 == len(coord_lst)):
            bv[i + 1] = 0
        else:
            currIndex += 1
        if (i % 1000000 == 0):
            print(i + 1)

    for i in range(coord_lst):
        bv[coord_lst[i] + 1] = 1

    return bv


## trial code
start_time = time.time()



print(bv[1003306:1003311])
print("--- %s seconds ---" % (time.time() - start_time))
    
# df.to_csv('update.csv', sep = "\t")




# def constructBitVectorArray(self, input):
#     #idea: create one bit vector for every chromosome, vector should be length 24, with 0 as a dummy element
#     #run through call helper
#     return #array of bit vector objects




# def writeFiles(self, input):
#     #write out the intersection based on bit vector??
#     #returns a text file

