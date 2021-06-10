from pcawg import *
import re


def readResults(result_file_path):
    with open(result_file_path) as f:
        lines = [line.rstrip() for line in f]
    results = {}
    for line in lines:
        if 'Record' in line:
            chr, pos = (get_chr_and_pos_from_record(line))
            if chr not in results:
                results[chr] = []
            results[chr].append(pos)
    return results

          

def get_chr_and_pos_from_record(record):
    delimiters = "(", "....", ")", ","
    regexPattern = '|'.join(map(re.escape, delimiters))
    split_record = re.split(regexPattern, record)
    chr = split_record[1]
    pos = split_record[2]
    return re.split("=", chr)[1], re.split("=", pos)[1]



def findMutationRate(result_file, coordinates_file):
    rates = {}




readResults("output_files/output_K562_distal_both_FDR_0.1.txt")


# asdf = essentialElementReadIn("data/K562_distal_both_FDR_0.1.txt")
# print(asdf)