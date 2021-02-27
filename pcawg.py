import pandas as pd
import os
import BitVector as bv
import sys

def essentialElementReadIn(file_path):
    #idea: read in file, split by chromosome, parse into 2D array of strings, where each subarray contains the coordinates of each chromosome

    df = pd.read_csv(file_path, sep = "\t", names = ['Chr_Num', 'Start', 'End'])

    #cleaning data
    df.Chr_Num = df.Chr_Num.str[3:] #getting rid of Chr tag

    df[['Start', 'End']] = df[["Start", 'End']].apply(pd.to_numeric) #changing start and end coords to ints

    chromosome_coord_map = {}

    for index, row in df.iterrows():
        chrom_num = row['Chr_Num']
    
        if chrom_num not in chromosome_coord_map.keys():
            chromosome_coord_map[chrom_num] = []
        chromosome_coord_map[chrom_num].extend([row['Start'], row['End']])


    return chromosome_coord_map

print(essentialElementReadIn("data/K562_distal_depleted_FDR_0.1.txt"))


#reading in data using pandas
# df = pd.read_csv("K562_distal_both_FDR_0.1.txt", sep = "\t", names = ['Chr_Num', 'Start', 'End'])

# #cleaning data
# df.Chr_Num = df.Chr_Num.str[3:] #getting rid of Chr tag

# df[['Start', 'End']] = df[["Start", 'End']].apply(pd.to_numeric) #changing start and end coords to ints

# chromosome_coord_map = {}

# for index, row in df.iterrows():
#     chrom_num = row['Chr_Num']
    
#     if chrom_num not in chromosome_coord_map.keys():
#         chromosome_coord_map[chrom_num] = []
#     chromosome_coord_map[chrom_num].extend([row['Start'], row['End']])


# print(chromosome_coord_map["X"])
# print(len(chromosome_coord_map["X"]))

    
# df.to_csv('update.csv', sep = "\t")




# def constructBitVectorArray(self, input):
#     #idea: create one bit vector for every chromosome, vector should be length 24, with 0 as a dummy element
#     #run through call helper
#     return #array of bit vector objects

# def constructBitVector(self, input):
#     #helper method
#     #idea: create bit vector initialized w/ 0's, run through files and initialize all values within coordinate range to 1's
#     #again, 0 is dummy value to minimize indexing confusion
#     return #a bit vector

# def writeFiles(self, input):
#     #write out the intersection based on bit vector??
#     #returns a text file

# #bad file read code
# #filename = 'K562_distal_depleted_FDR_0.1' # select file by type file name
# #with open(filename) as file_obj: # open file with open() function and make alias name from filename (file_obj)
#     #for line in file_obj: # using for looping to read entire content line-by-line
#         #print(line.rstrip())
        
# #potential write file code:
