from pcawg import *
import re


def readResults(result_file_path):
    """[summary]

    Args:
        result_file_path ([type]): [description]

    Returns:
        dict: a dictionary mapping each chromosome to a list of hits (pos) in that chromosome. Key (str): chromosme number; Value (int): position coordinate
    """
    with open(result_file_path) as f:
        lines = [line.rstrip() for line in f]
    results = {}
    for line in lines:
        if 'Record' in line:
            chr, pos = (get_chr_and_pos_from_record(line))
            if chr not in results:
                results[chr] = []
            results[chr].append(int(pos))
    return results

          

def get_chr_and_pos_from_record(record):
    """Given a Record entry (pyvcf format) from a results file, returns the chromosome number and the coordinate position

    Args:
        record: a record entry (pyvcf format)

    Returns:
        str, str: chromosome, position
    """
    delimiters = "(", "....", ")", ","
    regexPattern = '|'.join(map(re.escape, delimiters))
    split_record = re.split(regexPattern, record)
    chr = split_record[1]
    pos = split_record[2]
    return re.split("=", chr)[1], re.split("=", pos)[1]


def generate_list_of_coordinate_intervals(essentialElements):
    intervals = {}
    for chr in essentialElements:
        intervals[chr] = []
        for i in range(0,len(essentialElements[chr])-1, 2):
            interval = [essentialElements[chr][i], essentialElements[chr][i+1]]
            intervals[chr].append(interval)
    return intervals


def findMutationHits(results_file, coordinates_file):
    results = readResults(results_file)
    coordinates = essentialElementReadIn(coordinates_file)
    coordinate_intervals = generate_list_of_coordinate_intervals(coordinates)
    hits_by_interval = {}
    for chr in results:
        coordinates_on_chr = coordinate_intervals[chr]
        chr_results = {}
        for hit in results[chr]:
            for interval in coordinates_on_chr:
                if isInInterval(int(interval[0]), int(interval[1]), hit):
                    if str(interval) not in chr_results:
                        chr_results[str(interval)] = []
                    chr_results[str(interval)].append(hit)
        hits_by_interval[chr] = chr_results
    return hits_by_interval


def isInInterval(start, end, pos):
    return pos >= start and pos <= end


def findMutationRate(mutation_hits):
    rates = {}
    for chr in mutation_hits:
        chr_rates = {}
        for interval in mutation_hits[chr]:
            chr_rates[interval] = len(mutation_hits[chr][interval])
        rates[chr] = chr_rates
    return rates


def writeMutationRateOutput(file_name, mutation_rates):
    with open("mutation_rates_" + file_name, "w") as output_file:
        output_file.write('CHROMOSOME' + '\t' + 'INTERVAL' + '\t' + 'NUM_HITS' + '\t' + 'MUTATION_RATE' + '\n')
        for chr in mutation_rates:
            for interval in mutation_rates[chr]:
                output_file.write(chr + '\t' + interval + '\t' + str(mutation_rates[chr][interval]) + '\t' + str(mutation_rates[chr][interval]/getIntervalLength(interval)) + '\n')


def getIntervalLength(interval):
    delimiters = "[", "]", ","
    regexPattern = '|'.join(map(re.escape, delimiters))
    start_end = re.split(regexPattern, interval)
    return int(start_end[2]) - int(start_end[1]) + 1


for entry in os.scandir('data/essential_elements'):
    file_name = os.path.basename(entry)
    print(file_name)
    mutation_rate = findMutationRate(findMutationHits("output_files/output_" + file_name, "data/essential_elements/" + file_name))
    writeMutationRateOutput(file_name, mutation_rate)