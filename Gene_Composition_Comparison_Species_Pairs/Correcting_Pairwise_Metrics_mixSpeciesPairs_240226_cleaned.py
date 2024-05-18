# Correcting_Pairwise_Metrics_mixSpeciesPairs_240226.py

import pandas as pd
import numpy as np
from functools import reduce

####################
# Functions part I #
####################

# function 1: to remove the duplicate entries
# for BR - find the one that has the highest Jaccard value - likely the correct comparison with BR
# reason for this - BR will assign a gene to multiple clusters depending on its segments

def removeDuplicateEntries(subset):
    uniqueEntries = {}
    uniqueGene = []
    duplicateEntries = {}
    
    subset[["Jaccard"]] = subset[["Jaccard"]].apply(pd.to_numeric, errors='ignore')
    
    toRemove = subset[(subset["program2"] == "BR")]
    
    for index, row in toRemove.iterrows():
        #key1 = (row[2], row[3])
        
        speciesPair = row[1]
        anchorGene = row[2]
        OG2 = row[5]
        
        #print(speciesPair, anchorGene, OG2)
        
        if anchorGene not in uniqueGene:
            uniqueGene.append(anchorGene)

            geneSubset = toRemove[(toRemove["gene"] == anchorGene)]

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
                        duplicateEntries[key1] = set(anchorGene)
                    else:
                        duplicateEntries[key1].add(anchorGene)
                            
                        
    return(uniqueEntries, uniqueGene, duplicateEntries)

# function 2: check the size of the new dataframe without duplicate entries
def newMetricsDF(dictEntries):
    listOfDF = []
    for pairKey, metricRow in dictEntries.items():
        print(pairKey)
        columnNames = list(metrics_f.columns)
        progPairDF = pd.concat(metricRow)
        progPairDF_nodup = progPairDF.drop_duplicates(subset=progPairDF.columns.difference(['Unnamed: 0']))
        listOfDF.append(progPairDF_nodup)
        print(progPairDF.shape)
        print(progPairDF_nodup.shape)
    
    newMetrics = pd.concat(listOfDF)
    return(newMetrics)

# function 3: creating two dictionaries of dictionaries 
# dictionary of dictionaries
# dictionary 1 = program 1 and program 2 as keys, gene and all other info as values
# dictionary 2 = within each program 1 and program 2
# groups of arabidopsis genes (the union of p1 and p2) as keys, all other info as values (in a list)

def createDictionaries(subset):
    programPairsSame = {}
    programPairsDiff = {}
    countSameList = {}
    countDiffList = {}
    for index, row in subset.iterrows():
        spPair = row[1]
        prog2 = row[5]
        key1 = (spPair,prog2)
        
        anchor_gene = row[2]
        p1_num = row[-4]
        p2_num = row[-2]
        p1_anc = row[-3]
        p2_anc = row[-1]
        #file_to_write.write(str(p1_num) + "," + str(p2_num)+"\n")
        if p1_num == 0 and p2_num >= 1:
            anchor_union = p2_anc.split(";")
            anchor_union.sort()
            
            if key1 not in countDiffList:
                countDiffList[key1] = [anchor_gene]
            else:
                countDiffList[key1].append(anchor_gene)
                
            if key1 not in programPairsDiff:
                programPairsDiff[key1] = [set(anchor_union)]
                #file_to_write.write("1a: " + ATgene + "\n")
            else:
                count = 0
                for gene in anchor_union:
                    for geneSet in programPairsDiff[key1]:
                        if gene in geneSet:
                            count = count + 1
                            ind = programPairsDiff[key1].index(geneSet)
                            
                if count == 0:
                    programPairsDiff[key1].append(set(anchor_union))
                else: 
                    programPairsDiff[key1][ind].update(anchor_union)
                    #file_to_write.write("1c: "+ ATgene + "\n")
                
        elif p1_num >= 1 and p2_num == 0:
            anchor_union = p1_anc.split(",")
            anchor_union.sort()
            
            if key1 not in countDiffList:
                countDiffList[key1] = [anchor_gene]
            else:
                countDiffList[key1].append(anchor_gene)
                    
            if key1 not in programPairsDiff:
                programPairsDiff[key1] = [set(anchor_union)]
                #file_to_write.write("2a: "+ ATgene + "\n")
            else:
                count = 0
                for gene in anchor_union:
                    for geneSet in programPairsDiff[key1]:
                        if gene in geneSet:
                            count = count + 1
                            ind = programPairsDiff[key1].index(geneSet)
                            
                if count == 0:
                    programPairsDiff[key1].append(set(anchor_union))
                else: 
                    programPairsDiff[key1][ind].update(anchor_union)
                    #file_to_write.write("2c: "+ ATgene + "\n")
                
        elif (p1_num > 1) or (p2_num > 1):
            p1List = p1_anc.split(",") 
            p1List.sort()
            p2List = p2_anc.split(";")
            p2List.sort()
            if p1List == p2List:
                anchor_same = list(set(p1List + p2List))
                anchor_same.sort()
                
                if key1 not in countSameList: 
                    countSameList[key1] = [anchor_gene] 
                else:
                    countSameList[key1].append(anchor_gene)
                
                if key1 not in programPairsSame:
                    programPairsSame[key1] = [set(anchor_same)]
                    #file_to_write.write("3a: "+ ATgene + "\n")
                else:
                    programPairsSame[key1].append(set(anchor_same))
                    #file_to_write.write("3b: "+ ATgene+ "\n")
                        
            else:
                anchor_union = list(set(p1List + p2List))
                anchor_union.sort()
                
                if key1 not in countDiffList:
                    countDiffList[key1] = [anchor_gene]
                else:
                    countDiffList[key1].append(anchor_gene)
                    
                if key1 not in programPairsDiff:
                    programPairsDiff[key1] = [set(anchor_union)]
                    #file_to_write.write("4a: "+ ATgene+ "\n")
                else:
                    count = 0
                    for gene in anchor_union:
                        for geneSet in programPairsDiff[key1]:
                            if gene in geneSet:
                                count = count + 1
                                ind = programPairsDiff[key1].index(geneSet)
                    
                    if count == 0:
                        programPairsDiff[key1].append(set(anchor_union))
                    else: 
                        programPairsDiff[key1][ind].update(anchor_union)
                            #file_to_write.write("4c: "+ ATgene+ "\n")
    
    return(programPairsSame, programPairsDiff, countSameList, countDiffList)

# function 4: check the number of genes in the "keys" and number in the original list of genes
def checkNumGenes(dictionary):
    for key, item in dictionary.items():
        print(key)
        print(len(item))
        geneList =[]
        for geneSet in item:
            for gene3 in geneSet:
                geneList.append(gene3)
                
        print(len(geneList))
        print(len(set(geneList)))
        
# function 5: check if the extra genes/orthogroups for comparisons with BR are included in the keys
def checkGenesIncluded(dictionary1, dictionary2):
    extraGenes = {}
    for key, item in dictionary1.items():
        print(key)
        print(len(item))
        geneList =[]
        for geneSet in item:
            for gene in geneSet:
                geneList.append(gene)
        print(len(geneList))
        print(len(set(geneList)))
    
        genesNotInKey = set()
        for key2, item2 in dictionary2.items():            
            if key2 == key:
                print(key2)
                print(len(item))
                for gene2 in geneList:
                    if gene2 not in item2:
                        genesNotInKey.add(gene2)
        
        print(len(genesNotInKey))
        print(genesNotInKey)
        
        if len(genesNotInKey) > 0:
            if key not in extraGenes:
                extraGenes[key] = genesNotInKey
                
    return(extraGenes)

# function 6: for only the scenario where the AT genes are not exactly the same between the two orthogroups
# iteration/comparison among the "AT groups" to see if there are any duplicate genes between them
# stops only when only one final "AT group" to add to the list

def removeDup_newSet(testList):
    atgenes3 = []
    for at_tuple in testList:
        for tup in at_tuple:
            atgenes3.append(tup)

    #print(len(atgenes3))
    test_copy = testList.copy()

    newTupleList = testList

    while len(atgenes3) != len(set(atgenes3)):
        atgenes = []
        #print(atgenes)
        for i in range(len(newTupleList)):
            for j in range(len(newTupleList)-1):
                atgroup1 = newTupleList[i]
                atgroup2 = newTupleList[j+1]
                if atgroup1 != atgroup2:
                    if len(set(atgroup1).intersection(set(atgroup2))) > 0:
                        #print("at1")
                        #print(atgroup1)
                        #print("at2")
                        #print(atgroup2)
                        #print(atgroup1.intersection(atgroup2))
                        keyunion = list(set(atgroup1).union(set(atgroup2)))
                        keyunion.sort()
                        #print(keyunion)
                        if set(atgroup1) in test_copy:
                            test_copy.remove(set(atgroup1))
                        if set(atgroup2) in test_copy:
                            test_copy.remove(set(atgroup2)) 
                        if set(keyunion) not in test_copy:
                            test_copy.append(set(keyunion))
                        
                        atgenes.append(tuple(keyunion))
                     
                            
        atgenes_nodup = set(atgenes)
        #print(atgenes_nodup)
        
        atgenes3 = []
        for at_tuple in atgenes_nodup:
            for tup in at_tuple:
                atgenes3.append(tup)

        newTupleList = list(atgenes_nodup)

        #print("lengthTuple:")
        print(len(atgenes3))
        print(len(set(atgenes3)))
    
    print(len(testList))
    print(len(test_copy))
    return(test_copy)

# function 7: for only the scenario where the AT genes are not exactly the same between the two orthogroups
# removeDup_newSet function is found in this function
# updating the keys for each program pair

def updateKeys(dictionary):
    newDict = {}

    for programKey, keyList in dictionary.items():
        print(programKey)
        #print(keySet)
        geneList =[]
        for tup in keyList:
            for gene in tup:
                geneList.append(gene)
        print(len(geneList))
        print(len(set(geneList)))
        
        if len(geneList) == len(set(geneList)):
            newDict[programKey] = keyList
        else:
            newKeyList = removeDup_newSet(keyList)
            #print(newKeySet)
            newDict[programKey] = newKeyList
    return(newDict)

# function 8: for both the "same" and "diff" groups
# create a new dictionary with the new keys, especially the updated ones for the "diff" set

def createNewDictionary(dictionary):
    newDictionary = {}
    for programPair, atKeySet in dictionary.items():
        if programPair not in newDictionary:
            newDictionary[programPair] = {}
            for atKey in atKeySet:
                newATkey = list(atKey)
                newATkey.sort()
                newATkey = tuple(newATkey)
                if newATkey not in newDictionary[programPair]:
                    newDictionary[programPair][newATkey] = []
    return(newDictionary)

# function 9: add values to the dictionary of dictionaries
# values are individual rows
# if the AT gene is found in the AT key, which is the union of AT genes found between the two orthogroups
# metrics_cat_nodup, newProgramPairsDiff_f, extraBR_f
def addValuesToDictionaries_Diff(subset, dictionary, geneDict):    
    for index, row in subset.iterrows():
        spPair = row[1]
        prog2 = row[5]
        key1 = (spPair, prog2)
        gene = row[2]
        
        if key1 in geneDict:
            extraGenes = geneDict[key1]
            for anchor_key, anchor_items in dictionary[key1].items():
                #print(ATkey)
                if gene in anchor_key:
                    #print(gene)
                    if gene not in extraGenes:
                        dictionary[key1][anchor_key].append(row.tolist())
        else:
            for anchor_key, anchor_items in dictionary[key1].items():
                #print(ATkey)
                if gene in anchor_key:
                    dictionary[key1][anchor_key].append(row.tolist())
                
    return(dictionary) 

# function 10: finding the average metric values for each union of AT genes
# taking the average for all of the metrics
# however, hard to find the "average" for the number of AT genes in p1 and p2 and the AT gene compositions
# so taking the values from the first gene in the list

def newAvgList(programPairs):
    programKeys = programPairs.keys()
    newAvgValues = []
    programPairCount = []
    for pairKey in programKeys:
        print(pairKey)
        countList = []
        for anchorKey, anchorItems in programPairs[pairKey].items():
            geneOG = anchorKey
            #print(ATkey)
            program1 = anchorItems[0][3]
            program2 = anchorItems[0][5]
            intersectionAvg = sum(float(x[7]) for x in anchorItems)/len(anchorItems)
            unionAvg = sum(float(x[8]) for x in anchorItems)/len(anchorItems)
            JIavg = sum(float(x[9]) for x in anchorItems)/len(anchorItems)
            p1_num = anchorItems[0][10]
            p1_anchor_list = anchorItems[0][11]
            p2_num = anchorItems[0][12]
            p2_anchor_list = anchorItems[0][13]
            #print([geneOG, program1, program2, RSavg, ARSavg, intersectionAvg, unionAvg,JIavg,p1_num,p1_ATlist,p2_num,p2_ATlist])
            newAvgValues.append([pairKey[0],geneOG, program1, program2, intersectionAvg, unionAvg,JIavg,p1_num,p1_anchor_list,p2_num,p2_anchor_list])
            countList.append(len(anchorKey))
            
        finalCount = sum(countList) 
        print(finalCount)
        programPairCount.append([pairKey, finalCount])
        
    return(newAvgValues, programPairCount)


#####################
# Functions Part II #
#####################

# function 1: get the initial summary numbers from the original metrics table
def getSummaryNumbers(metrics):
    programPairs= {}
    for index, row in metrics.iterrows():
        spPair = row[1]
        prog2 = row[5]
        key1 = (spPair,prog2)
        
        if key1 not in programPairs:
            programPairs[key1] = []
            programPairs[key1].append(row.tolist())
        else:
            programPairs[key1].append(row.tolist())
            
    sumList = []
    for pairKey, metricRow in programPairs.items():
        print(pairKey)
        columnNames = list(metrics.columns)
        progPairDF = pd.DataFrame(metricRow, columns = columnNames)
        print(progPairDF.shape)
        #print(progPairDF)
        sub1 = progPairDF[(progPairDF["p1_num_anc"] == 1) & (progPairDF["p2_num_anc"] == 1)]
        sub1 = sub1.drop_duplicates(subset=sub1.columns.difference(['Unnamed: 0']))
        print(sub1.shape)
        sub2 = progPairDF[(progPairDF["p1_num_anc"] > 1) | (progPairDF["p2_num_anc"] > 1)]
        sub2 = sub2.drop_duplicates(subset=sub2.columns.difference(['Unnamed: 0']))
        print(sub2.shape)
        sub3 = progPairDF[(progPairDF["p1_num_anc"] == 0) & (progPairDF["p2_num_anc"] == 1)]
        sub3 = sub3.drop_duplicates(subset=sub3.columns.difference(['Unnamed: 0']))
        #print(sub3)
        print(sub3.shape)
        sub4 = progPairDF[(progPairDF["p1_num_anc"] == 1) & (progPairDF["p2_num_anc"] == 0)]
        sub4 = sub4.drop_duplicates(subset=sub4.columns.difference(['Unnamed: 0']))
        #print(sub4)
        print(sub4.shape)
        
        fulldiff_anc = sub2.shape[0] + sub3.shape[0] + sub4.shape[0]
    
        sumList.append([pairKey, progPairDF.shape[0], sub1.shape[0], fulldiff_anc])
    
    summaryDF = pd.DataFrame(sumList, columns = ["spPair_program2","All","one_anc","diff_anc"])
    summaryDF["Total"] = summaryDF["one_anc"] + summaryDF["diff_anc"]
    summaryDF["PropOne"] = summaryDF["one_anc"]/summaryDF["Total"]
    
    return(summaryDF)

# function 2: counting the number of orthogroups for comparison from the updated list

def countrowsNewDF(dataframe):
    programPairs= {}
    for index, row in dataframe.iterrows():
        spPair = row[0]
        prog2 = row[3]
        key1 = (spPair,prog2)
        
        if key1 not in programPairs:
            programPairs[key1] = []
            programPairs[key1].append(row.tolist())
        else:
            programPairs[key1].append(row.tolist())
    
    countList = []
    for pairKey, metricRow in programPairs.items():
        print(pairKey)
        columnNames = list(dataframe.columns)
        progPairDF = pd.DataFrame(metricRow, columns = columnNames)
        progPairDF_nodup = progPairDF.drop_duplicates(subset=progPairDF.columns.difference(['Unnamed: 0']))
        print(progPairDF.shape)
        print(progPairDF_nodup.shape)
        countList.append([pairKey,progPairDF_nodup.shape[0]])
    
    return(countList)

# function 3: combine all the outputs into a summary dataframe

def combineDF(metricsDF, ori_cts1_d, newDiffDF):
    initialSummary = getSummaryNumbers(metricsDF)
    
    oriCountDiff_DF = pd.DataFrame(ori_cts1_d, columns=["spPair_program2","oriATgenes_diff"])
    
    list_diff = countrowsNewDF(newDiffDF)
    
    newCountDiff_DF = pd.DataFrame(list_diff, columns=["spPair_program2","combAvg_diff"])
    
    listOfDF = [initialSummary, oriCountDiff_DF, newCountDiff_DF]
    df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['spPair_program2']), listOfDF)
    
    df_merged["newTotalOG"] = df_merged["one_anc"] + df_merged["combAvg_diff"]
    df_merged["oriTotal"] = df_merged["oriATgenes_diff"] 
    
    return(df_merged)

#########################################################################
# Functions Part 3 - get summary JI for each species pair and algorithm #
#########################################################################

def calcJaccardList(mergeDF):
    spPair_prog_List = []
    sumList = []
    
    for index, row in mergeDF.iterrows():
        spPair = row[0]
        prog2 = row[4]
        
        if [spPair, prog2] not in spPair_prog_List:
            spPair_prog_List.append([spPair, prog2])
            
            subset = mergeDF[(mergeDF["speciesPair"] == spPair) & (mergeDF["program2"] == prog2)]
        
            
            countOne = subset["Jaccard"].value_counts()[1]
            totalOG = len(subset)
            propOne = countOne / totalOG
            meanJI = subset["Jaccard"].mean()
            sdJI = subset["Jaccard"].std()
            seJI = subset["Jaccard"].sem()
           
            print(prog2,spPair,countOne,totalOG,propOne,meanJI,sdJI,seJI)
            sumList.append([prog2,spPair,countOne,totalOG,propOne,meanJI,sdJI,seJI])
        
        columnNames = ["algorithm","SpeciesPair","numberOfExactOG","totalOG","proportionOfExactOG",
                       "meanJI","sdJI","seJI"]
        newJISumDF = pd.DataFrame(sumList, columns = columnNames)
    
    return(newJISumDF)  

#########################
# START HERE with files #
#########################

# files
main_path = r"/path_to_files"
full_file = main_path + "/pairs/output/mixed_speciesPairsJIvalues_toCorrect_240226.csv"

metrics_f = pd.read_csv(full_file, dtype={"OG2": str})
metrics_f.replace('na', np.nan, inplace=True)
metrics_f.shape

print("remove duplicate entries")
dictEntries_f, uniqueList_f, duplicateDict_f = removeDuplicateEntries(metrics_f)
print("duplicate entries")
for pairKey, metricRow in duplicateDict_f.items():
    print(pairKey)
    print(metricRow)
    print(len(metricRow))

# create updated metrics table with no duplicate gene orthogroup entries
print("create new metrics table with no duplicate entries")    
toKeep = metrics_f[(metrics_f["program2"] != "BR")]
metrics_full_BR_nodup = newMetricsDF(dictEntries_f)
metrics_cat = pd.concat([toKeep, metrics_full_BR_nodup])
metrics_cat_nodup = metrics_cat.drop_duplicates()
metrics_cat_nodup.shape

# create dictionary of keys 
print("create Dictionary")
sameATnum_f, diffATnum_f, countSame_f, countDiff_f = createDictionaries(metrics_cat_nodup)

# assessing duplicate entries in orthogroup comparisons where the number of anchor genes are not the same
for key, item in countDiff_f.items():
    print(key)
    print(len(item))
    print(len(set(item)))
    uniqueList = []
    dupList = []
    for gene in item:
        if gene not in uniqueList:
            uniqueList.append(gene)
        else: 
            dupList.append(gene)
    print(dupList)

# check gene numbers
checkNumGenes(diffATnum_f)

# check if any extra BR genes, which ones are they. 
# Will make sure not to include them when adding values to the new dictionary
extraBR_f = checkGenesIncluded(diffATnum_f, countDiff_f)

# for gene/orthogroup comparisons with different number of anchor genes between them
print("iterative method to get equal genes in dictionary keys of anchor genes")
noDupGeneKey = updateKeys(diffATnum_f)
print("check number of genes again")
checkNumGenes(noDupGeneKey)
print("create new dictionary of dictionaries with new AT keys")
newProgramPairsDiff_f = createNewDictionary(noDupGeneKey)
print("add values")
diffATvalue_f = addValuesToDictionaries_Diff(metrics_cat_nodup, newProgramPairsDiff_f, extraBR_f)
print("calc new avg")
diffOGList_f, oriATnum_diff_f = newAvgList(diffATvalue_f)

#create dataframes for both datasets
columnNames = ["speciesPair","gene","program1","program2","intersection","union","Jaccard","p1_num_anc",
                   "p1_anc_genes","p2_num_anc","p2_anc_genes"]
diffOGdf_f = pd.DataFrame(diffOGList_f, columns=columnNames)

#############################################
# call function to create summary dataframe #
#############################################
full_sum = combineDF(metrics_cat_nodup, oriATnum_diff_f, diffOGdf_f)
summaryOutput_f = main_path + "/pairs/output/summary_mixedPairs_240227.csv"
full_sum.to_csv(summaryOutput_f, index = False)

# Creating the final updated metrics table for making heatmaps
oneSubset_f = metrics_cat_nodup[(metrics_cat_nodup["p1_num_anc"]== 1) & (metrics_cat_nodup["p2_num_anc"] == 1)]
oneSubset_f = oneSubset_f.drop(["Unnamed: 0"], axis=1)
print(oneSubset_f.shape)
print("concatenate dataframes")
metrics_dfs_f = [oneSubset_f, diffOGdf_f]
metrics_merged_f = pd.concat(metrics_dfs_f)
print(metrics_merged_f.shape)
print("drop duplicates")
nodup_merged_f = metrics_merged_f.drop_duplicates()
print(nodup_merged_f.shape)
newMetricsNoDup_out_f = main_path + "/pairs/output/spPairs_mixed_newMetrics_noDup_240227.csv"
nodup_merged_f.to_csv(newMetricsNoDup_out_f, index=False)

########################################################
# get summary JI for each species pair and algorithm 2 #
########################################################
mixJIsummary = calcJaccardList(nodup_merged_f)
mixJIsummary_out = main_path + "/pairs/output/spPairs_mixed_JIsummary_240227.csv"
mixJIsummary.to_csv(mixJIsummary_out, index=False)