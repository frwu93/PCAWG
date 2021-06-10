from pcawg import *

## trial code
start_time = time.time()

chromosome_coord_map_essential_both = essentialElementReadIn("data/K562_distal_both_FDR_0.1.txt")
chromosome_coord_map_essential_depleted = essentialElementReadIn("data/K562_distal_depleted_FDR_0.1.txt")
chromosome_coord_map_essential_enriched = essentialElementReadIn("data/K562_distal_enriched_FDR_0.1.txt")
chromosome_coord_map_essential_nonsig = essentialElementReadIn("data/K562_distal_nonsig_FDR_0.1.txt")

chromosome_coord_map_0g = essentialElementReadIn("data/distalDHS_0gRNAs_FDR0.1.txt")
chromosome_coord_map_1g = essentialElementReadIn("data/distalDHS_1gRNAs_FDR0.1.txt")

chromosome_length_map = chromosomeLengthReadIn("data/ChromosomeLengths.txt")

bvMapEssential_both = constructBitVectorMap(chromosome_coord_map_essential_both, chromosome_length_map)
bvMapEssential_depleted = constructBitVectorMap(chromosome_coord_map_essential_depleted, chromosome_length_map)
bvMapEssential_enriched = constructBitVectorMap(chromosome_coord_map_essential_enriched, chromosome_length_map)
bvMapEssential_nonsig = constructBitVectorMap(chromosome_coord_map_essential_nonsig, chromosome_length_map)

bvMap0g = constructBitVectorMap(chromosome_coord_map_0g, chromosome_length_map)
bvMap1g = constructBitVectorMap(chromosome_coord_map_1g, chromosome_length_map)

map_of_maps = {'K562_distal_both_FDR_0.1.txt':bvMap0g, 'K562_distal_depleted_FDR_0.1.txt':bvMap1g,
'K562_distal_enriched_FDR_0.1.txt':bvMapEssential_enriched, 'K562_distal_nonsig_FDR_0.1.txt': bvMapEssential_nonsig}

no_result_maps = {'distalDHS_0gRNAs_FDR0.1.txt':bvMap0g, 'distalDHS_1gRNAs_FDR0.1.txt':bvMap1g}

results = compare("data/final_consensus_decompressed/snv_mnv", no_result_maps)
for bvMap in results:
    with open("output_" + bvMap, "w") as output_file:
        for entry in results[bvMap]:
            output_file.write(os.path.basename(entry) + '\n')
            for line in results[bvMap][entry]:
                info = ''.join(key + str(val) + ';' for key, val in line[1].items())
                output_file.write(line[0] + '....' + info + '\n')



# print(len(results))
print("--- %s seconds ---" % (time.time() - start_time))
    
# df.to_csv('update.csv', sep = "\t")
