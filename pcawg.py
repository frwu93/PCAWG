import pandas as pd
import os
from BitVector import BitVector
import sys
import time 


def essentialElementReadIn(file_path):
    """Reads in essential element coordinates from TSV file and returns a map containing the information. Maps chromosome number
    to a list of that chromosome's essential element coordinates.

    Args:
        file_path (str): file path of desired file

    Returns:
        dict: Dictionary mapping each chromosome to its essential element coordinates. Key: (str) chromosome number. 
        Value: (list) list of essential element coordinates in format of [start1, end1, start2, end2, ... , start[n], end[n]]
    """
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
    """Reads in chromosome lengths from CSV file and returns a map containing the information. Maps chromosome number to 
    its corresponding length.

    Args:
        file_path (str): file path of desired file

    Returns:
        dict: Dictionary mapping each chromosome to its corresponding length. Key: (str) chromosome number. Value: (int) length of chromosome in bp
    """
    df = pd.read_csv(file_path, names = ['Chr_Num', 'Length'])
    df['Length'] = pd.to_numeric(df['Length']) 

    chromosome_length_map = {}

    for index, row in df.iterrows():
        chrom_num = row['Chr_Num']
        chromosome_length_map[chrom_num] = row['Length']

    return chromosome_length_map


def constructBitVector(chromosome_coord_map, chromosome_length_map, chr_num):
    """Constructs a bit vector for an individual chromosome.

    Args:
        chromosome_coord_map (dict): map of chromosomes and desired coordinates. Key: (str) chromosome number. Value: (list) Desired coordinates 
        chromosome_length_map (dict): map of chromosomes to their lengths. Key: (str) chromosome number. Value: (int) Length of chromosome. 
        chr_num (str): chromosome number 

    Returns:
        BitVector: a bit vector representing regions in chromosome that exhibited unusual growth activity.
        0 for coordinates not associated with cell growth, 1 for coordinates associated with cell growth
    """  
    currIndex = 0

    bv  = BitVector(size = chromosome_length_map[chr_num] + 1) # initializes a bit vector of 0's of size of chromosome + 1

    coord_lst = chromosome_coord_map[chr_num]

    for i in range(chromosome_length_map[chr_num]):
        if (i < coord_lst[currIndex] and currIndex % 2 == 0): # previous end coordinate < index <  current start coordinate
            bv[i] = 0
        elif(i < coord_lst[currIndex] and currIndex % 2 == 1): # current start coordinate < index <  current end coordinate
            bv[i] = 0
            bv[i] = 1
        elif(currIndex + 1 == len(coord_lst)): # very last end coordinate < index < end of chromosome
            bv[i] = 0
        else: # index >= current end coordinate --> we have entered a new region
            currIndex += 1
        if (i % 1000000 == 0): # testing purposes
            print(i)

    for i in range(len(coord_lst)): # go back and initialize all boundries to 1
        bv[coord_lst[i]] = 1

    return bv


def constructBitVectorMap(chromosome_coord_map, chromosome_length_map):
    """Constructs and returns a map of chromosomes and their corresponding bit vectors. Each bit vector represents a chromosome's
    essential element coordinates. Key: (str) chromosome number. Value: (BitVector) bit vector representation of coordinates for that chromosome

    Args:
       chromosome_coord_map (dict): map of chromosomes and desired coordinates. Key: (str) chromosome number. Value: (list) Desired coordinates 
       chromosome_length_map (dict): map of chromosomes to their lengths. Key: (str) chromosome number. Value: (int) Length of chromosome. 

    Returns:
        dict: map of chromosomes to their corresponding bit vector
    """

    bitVectorMap = {}

    for chr_num in chromosome_coord_map.keys():
        bv = constructBitVector(chromosome_coord_map, chromosome_length_map, chr_num)
        bitVectorMap[chr_num] = bv
   
    return bitVectorMap

## trial code
start_time = time.time()

chromosome_coord_map = essentialElementReadIn("data/K562_distal_both_FDR_0.1.txt")

chromosome_length_map = chromosomeLengthReadIn("data/ChromosomeLengths.txt")

bv = constructBitVector(chromosome_coord_map, chromosome_length_map, "22")


print(bv[19784619:19784624])
print(bv[19785041:19785046])

print("--- %s seconds ---" % (time.time() - start_time))
    
# df.to_csv('update.csv', sep = "\t")









# def writeFiles(self, input):
#     #write out the intersection based on bit vector??
#     #returns a text file

