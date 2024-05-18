# CompareOrthogroupComposition_SpeciesPairs_dip_240222.py
# new code for calculating jaccard index for species pairs compared to baseline

import pandas as pd
import numpy as np
from sklearn.metrics.cluster import rand_score
from sklearn.metrics.cluster import adjusted_rand_score
import itertools
import re
import os

################################################
# function 1: read in files                    #        
#             sort genes in order within cells #
#             change the column headers        #
################################################
def sortGenes(OGfile, program):
    
    orthogroup = pd.read_table(OGfile)
    print("orthogroup shape:", orthogroup.shape)
    OGsort = orthogroup.fillna("*")
    
    if program == "orthofinder":
        
        for row, value in OGsort.iterrows():
            index = OGsort.shape[1]

            i=0
            while i < index:
                items = value[i]
                #print(orthogroup.iloc[row][i])
                if type(items) != int:
                    geneList = items.split(", ")
                    #print(geneList)
                    if len(geneList) > 1:
                        #print("Old:", OGsort.iat[row,i])
                        geneList.sort()
                        #print(geneList)
                        sortedGenes = ','.join(geneList)
                        #print(sortedGenes)
                        OGsort.iat[row,i] = sortedGenes
                        #print("New:", OGsort.iat[row,i])
                        i += 1
                    else:
                        #print("Length 1:", OGsort.iloc[row][i])
                        i += 1
                else:
                    i += 1 
            
    elif program == "sonicparanoid":
        
        OGsort.columns = [col_name.split(".")[0] for col_name in OGsort.columns]
        
        for row, value in OGsort.iterrows():
            index = OGsort.shape[1]
            
            i=0
            while i < index:
                items = value[i]
                #print(orthogroup.iloc[row][i])
                if type(items) != int:
                    geneList = items.split(",")
                    #print(geneList)
                    if len(geneList) > 1:
                        #print("Old:", OGsort.iat[row,i])
                        geneList.sort()
                        #print(geneList)
                        sortedGenes = ','.join(geneList)
                        #print(sortedGenes)
                        OGsort.iat[row,i] = sortedGenes
                        #print("New:", OGsort.iat[row,i])
                        i += 1
                    else:
                        #print("Length 1:", OGsort.iloc[row][i])
                        i += 1
                else:
                    i += 1
                    
    elif program == "broccoli":
        OGsort.columns = [col_name.split(".")[0] for col_name in OGsort.columns]
        OGsort = OGsort.rename(columns={"#OG_name": "OG_name"})
        
        for row, value in OGsort.iterrows():
            index = OGsort.shape[1]

            i=0
            while i < index:
                items = value[i]
                #print(orthogroup.iloc[row][i])
                if type(items) != int:
                    geneList = items.split(" ")
                    #print(geneList)
                    if len(geneList) > 1:
                        #print("Old:", OGsort.iat[row,i])
                        geneList.sort()
                        #print(geneList)
                        sortedGenes = ','.join(geneList)
                        #print(sortedGenes)
                        OGsort.iat[row,i] = sortedGenes
                        #print("New:", OGsort.iat[row,i])
                        i += 1
                    else:
                        #print("Length 1:", OGsort.iloc[row][i])
                        i += 1
                else:
                    i += 1
                    
    elif program == "orthnet":
        for row, value in OGsort.iterrows():
            index = OGsort.shape[1]

            i=0
            while i < index:
                items = value[i]
                #print(orthogroup.iloc[row][i])
                if type(items) != int:
                    geneList = items.split(", ")
                    #print(geneList)
                    if len(geneList) > 1:
                        #print("Old:", OGsort.iat[row,i])
                        geneList.sort()
                        #print(geneList)
                        sortedGenes = ','.join(geneList)
                        #print(sortedGenes)
                        OGsort.iat[row,i] = sortedGenes
                        #print("New:", OGsort.iat[row,i])
                        i += 1
                    else:
                        #print("Length 1:", OGsort.iloc[row][i])
                        i += 1
                else:
                    i += 1     
    return(OGsort)

######################################
# function 2: assign cluster to gene #
######################################
def clusterGeneList(OGdf, program):
    
    # create first dictionary with cluster as key, and list of all the genes in the cluster as the items
    clusterList = []
    for index, row in OGdf.iterrows():
        #print(row)
        cluster = row[0] 
        #print(cluster)
        
        valueList = list(row)
        naCount = valueList.count("*")
        if naCount < 4:
            for i in range(len(row)-1):
                geneList = row[i+1]
                geneList = geneList.split(",")
                for gene in geneList:
                    if gene != "*":
                        clusterList.append([gene, cluster])
    
    clusterDF = pd.DataFrame(clusterList, columns=["gene", program])
    return(clusterDF)
    
#####################################################
# function 3: assign cluster to gene - for baseline #
#####################################################
def clusterGeneList_baseline(OGdf, program):
    
    newOF = OGdf.fillna("*")
    # create first dictionary with cluster as key, and list of all the genes in the cluster as the items
    clusterList = []
    for index, row in newOF.iterrows():
        #print(row)
        cluster = row[0] 
        #print(cluster)
        
        valueList = list(row)
        naCount = valueList.count("*")
        if naCount < 1:
            for i in range(len(row)-1):
                geneList = row[i+1]
                geneList = geneList.split(",")
                for gene in geneList:
                    if gene != "*":
                        clusterList.append([gene, cluster])
    
    clusterDF = pd.DataFrame(clusterList, columns=["gene", program])
    return(clusterDF)

######################################################################
# function 4: wrapper function for first 2 functions                 #
#             returns a dataframe with genes linked with the cluster #
######################################################################
def methodsDF(filePath, OGprogram, columnOG, sp1, sp2):
    sortDF = sortGenes(filePath, OGprogram)
    
    if OGprogram == "orthofinder":
        spPairDF = sortDF[["HOG", sp1, sp2]]
    elif OGprogram == "sonicparanoid":
        spPairDF = sortDF[["group_id", sp1, sp2]]
    elif OGprogram == "broccoli":
        spPairDF = sortDF[["OG_name",sp1, sp2]]
    elif OGprogram == "orthnet":
        spPairDF = sortDF[["groupNumber", sp1, sp2]]
    
    cluster_df = clusterGeneList(spPairDF, columnOG)
    return(cluster_df)

############################################################################
# function 5: to determine most inclusive OG composition for each scenario #
############################################################################
def maxSubsetGroup(subsetList, OG_df): 
    programMax = max(subsetList, key=lambda x: x[2])[0]
    clusterMax = max(subsetList, key=lambda x: x[2])[1]
    subset = OG_df.loc[OG_df[programMax]==clusterMax]
    subset = subset.fillna("*")
    return(subset)

#######################################################
# function 6: calculate the Jaccard Index for each OG # 
#######################################################
def calcJaccardList(mergeDF, columnOG, geneStart, speciesPair):
    metricsList = []
    clusterList = set()
    for index, row in mergeDF.iterrows():
        gene = row[0]
        if gene.startswith(geneStart):
            print(gene)

            p1 = "OFbase"
            p2 = columnOG
            c1 = row[1] #cluster1, OF baseline
            c2 = row[2] #cluster2, other program "columnOG"

            # calculate the jaccard proportion similarity - find the same gene to anchor
            
            p1_only = set(mergeDF.loc[mergeDF[p1]==c1]["gene"].tolist())
            p2_only = set(mergeDF.loc[mergeDF[p2]==c2]["gene"].tolist())
            
            # determine how many of anchor (anc) species genes are in each orthogroup
            p1_anc = []
            for p1_gene in p1_only:
                if p1_gene.startswith(geneStart):
                    p1_anc.append(p1_gene)
            p1_anc_len = len(p1_anc)
            p1_anc_join = ",".join(p1_anc)
                
            p2_anc = []
            for p2_gene in p2_only:
                if p2_gene.startswith(geneStart):
                    p2_anc.append(p2_gene)
            p2_anc_len = len(p2_anc)
            p2_anc_join = ";".join(p2_anc)
            
            # calculate JI
            intersect = len(p1_only & p2_only)
            un = len(p1_only) + len(p2_only) - intersect
            print(p1_only, p2_only)
            if intersect != 0:
                jaccard = float((intersect) / un)
            else:
                jaccard = 0
            
            print(jaccard)
            metricsList.append([speciesPair,gene,p1, c1, p2, c2, intersect,un,jaccard,p1_anc_len,p1_anc_join,p2_anc_len,p2_anc_join])
            
    #countOne = metricsList.count(1)
    #totalOG = len(metricsList)
    #propOne = countOne / totalOG
    #meanJI = sum(metricsList) / totalOG
    #sdJI = np.std(metricsList)
    #seJI = sdJI / np.sqrt(np.size(metricsList))
    #print(speciesPair,countOne,totalOG,propOne,meanJI,sdJI,seJI)
    
    metricColumns = ["speciesPair","gene", "program1", "OG1", "program2", "OG2", "intersection", "union", "Jaccard","p1_num_anc","p1_anc_genes","p2_num_anc","p2_anc_genes"]
    metricsDF = pd.DataFrame(metricsList, columns = metricColumns)
    
    return(metricsDF)

######################################################################################
# function 7: master wrapper for each method to test each species pair - diploid set #
######################################################################################

def masterSpeciesPairs(filePath, OGprogram, columnOG):
    speciesList = ['Ath', 'Chi', 'Cru', 'Tar', 'Aar']
    speciesPairMasterList = []
    pairs = itertools.combinations(speciesList, 2)
    for item in pairs:
        species1 = item[0]
        species2 = item[1]
        spPair = species1 + species2
        
        pairFile = r"/u/home/i/irenelia/orthology/pairs/"+species1 + species2+"/Orthogroups/Orthogroups.tsv"
        pair_df = pd.read_table(pairFile)
        speciesPair_df = clusterGeneList_baseline(pair_df, "OFbase")
        diploid_df = methodsDF(filePath, OGprogram, columnOG, species1, species2)
        mergeDF = speciesPair_df.merge(diploid_df, on="gene", how="outer")
        
        if species1 == "Ath":
            JI_DF = calcJaccardList(mergeDF, columnOG, "AT", spPair)
            speciesPairMasterList.append(JI_DF)
        elif species1 == "Chi":
            JI_DF = calcJaccardList(mergeDF, columnOG, "CARHR", spPair)
            speciesPairMasterList.append(JI_DF)
        elif species1 == "Cru":
            JI_DF = calcJaccardList(mergeDF, columnOG, "Carub", spPair)
            speciesPairMasterList.append(JI_DF)
        elif species1 == "Tar":
            JI_DF = calcJaccardList(mergeDF, columnOG, "gene", spPair)
            speciesPairMasterList.append(JI_DF)
    
    #columnNames = ["speciesPair","gene", "program1", "OG1", "program2", "OG2", "intersection", "union", "Jaccard","p1_num_anc","p1_anc_genes","p2_num_anc","p2_anc_genes"]
    sum_df = pd.concat(speciesPairMasterList)
    sum_df_nodup = sum_df.drop_duplicates()
    
    return(sum_df_nodup)  

################
################
## START HERE ##
################
################

#working directory
main_path = "/path_to_files"

print("Reading in files")
#paths to orthology output files
broccoli = main_path + "/dip_BR-table_OGs_protein_names.txt"
orthofinder_b = main_path + "/dip_OF_b_N0.tsv"
orthofinder_d = main_path + "/dip_OF_d_N0.tsv"
orthofinder_m = main_path + "/dip_OF_m_N0.tsv"
orthnet = main_path + "/OrthNet_groups_diploid_230927.txt"
sonicparanoid_d = main_path + "/fdip_SP_d-flat.ortholog_groups.tsv"
sonicparanoid_m = main_path + "/dip_SP_m-flat.ortholog_groups.tsv"

# call wrapper function
BR_df = masterSpeciesPairs(broccoli, "broccoli", "BR")
OFb_df = masterSpeciesPairs(orthofinder_b, "orthofinder", "OFb")
OFd_df = masterSpeciesPairs(orthofinder_d, "orthofinder", "OFd")
OFm_df = masterSpeciesPairs(orthofinder_m, "orthofinder", "OFm")
ON_df= masterSpeciesPairs(orthnet, "orthnet", "ON")
SPd_df = masterSpeciesPairs(sonicparanoid_d, "sonicparanoid", "SPd")
SPm_df = masterSpeciesPairs(sonicparanoid_m, "sonicparanoid", "SPm")

# combine all the summary dataframes, write output
methodsDF = [BR_df, OFb_df, OFd_df, OFm_df, ON_df, SPd_df, SPm_df]
methods_concat = pd.concat(methodsDF)
methods_concat_nodup = methods_concat.drop_duplicates()

summary_path = "/path_to_output/diploid_speciesPairsJIvalues_toCorrect_240222.csv"
methods_concat.to_csv(summary_path, sep = ",", index="False")  


