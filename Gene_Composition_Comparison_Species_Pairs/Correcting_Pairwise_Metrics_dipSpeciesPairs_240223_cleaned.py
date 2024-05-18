# Correcting_Pairwise_Metrics_dipSpeciesPairs_240223.py

# Part I: 
import pandas as pd
import numpy as np
from functools import reduce

# function 1: to remove the duplicate entries
# for BR - find the one that has the highest Jaccard value - likely the correct comparison with BR
# reason for this - BR will assign a gene to multiple clusters depending on its segments

def removeDuplicateEntries(subset):
    uniqueEntries = {}
    uniqueGene = []
    duplicateEntries = {}
    
    subset[["Jaccard"]] = subset[["Jaccard"]].apply(pd.to_numeric, errors='ignore')
    
    for index, row in subset.iterrows():
        #key1 = (row[2], row[3])
        
        speciesPair = row[1]
        anchorGene = row[2]
        OG2 = row[5]
        
        if OG2 == "BR":
            if anchorGene not in uniqueGene:
                uniqueGene.append(anchorGene)

                geneSubset = subset[(subset["gene"] == anchorGene)]

                for index2, row2 in geneSubset.iterrows():
                    spPair2 = row2[1]
                    
                    key1 = spPair2 #SpeciesPair
                
                    keySubset = geneSubset[(geneSubset["speciesPair"] == spPair2)]
                
                    if keySubset.shape[0] == 1:
                        if key1 not in uniqueEntries:
                            uniqueEntries[key1] = [keySubset]
                        else:
                            uniqueEntries[key1].append(keySubset)

                    else:
                        #print(keySubset)
                        #keySubset[["Jaccard"]] = keySubset[["Jaccard"]].apply(pd.to_numeric, errors='ignore')
                        maxJaccard = keySubset.nlargest(1,['Jaccard'])
                        #print(maxJaccard)
                        if key1 not in uniqueEntries:
                            uniqueEntries[key1] = [maxJaccard]
                        else:
                            uniqueEntries[key1].append(maxJaccard)

                        if key1 not in duplicateEntries:
                            duplicateEntries[key1] = [anchorGene]
                        else:
                            duplicateEntries[key1].append(anchorGene)
                            
        else:
            if key1 not in uniqueEntries:
                uniqueEntries[key1] = [row]
            else:
                uniqueEntries[key1].append(row)
                        
    return(uniqueEntries, uniqueGene, duplicateEntries)

    
# function 2: check the size of the new dataframe without duplicate entries
# concatenate rows and get new dataframe
def newMetricsDF(dictEntries, metricsDF):
    listOfDF = []
    for spPairKey, metricRow in dictEntries.items():
        print(spPairKey)
        columnNames = list(metricsDF.columns)
        progPairDF = pd.concat(metricRow)
        progPairDF_nodup = progPairDF.drop_duplicates(subset=progPairDF.columns.difference(['Unnamed: 0']))
        listOfDF.append(progPairDF_nodup)
        print(progPairDF.shape)
        print(progPairDF_nodup.shape)
    
    newMetrics = pd.concat(listOfDF)
    return(newMetrics)  

# files
main_path = r"/path_to_files"
#test_file = main_path + "/comparison/240209_full/test_metrics.csv"
dip_file = main_path + "/pairs/output/diploid_speciesPairsJIvalues_toCorrect_240222.csv"
#full_file = main_path + "/comparison/240209_full/full_pairwiseMetrics_240209.csv"

metrics_d = pd.read_csv(dip_file, dtype={"OG2": str})
#metrics_d.replace('na', np.nan, inplace=True)
print(metrics_d.shape)

# create files for updated metrics file
print("remove duplicate entries")
dictEntries_d, uniqueList_d, duplicateDict_d = removeDuplicateEntries(metrics_d)

# double check which ones have duplicate entries
print("double check which entries have duplicate genes")
for pairKey, metricRow in duplicateDict_d.items():
    print(pairKey)
    print(metricRow)

# create new metrics dataframe
metrics_dip_nodup = newMetricsDF(dictEntries_d, metrics_d)
print(metrics_dip_nodup.shape)

no_dup_metrics_output = main_path + "/pairs/output/diploid_speciesPairsJIvalues_noDuplicates_240223.csv"
metrics_dip_nodup.to_csv(no_dup_metrics_output, index=False)

