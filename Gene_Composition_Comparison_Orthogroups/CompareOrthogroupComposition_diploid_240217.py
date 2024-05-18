# CompareOrthogroupComposition_diploid_FINAL.py
# new code for calculating metrics (rand score, adjusted rand score, jaccard index)

import pandas as pd
import numpy as np
from sklearn.metrics.cluster import rand_score
from sklearn.metrics.cluster import adjusted_rand_score
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
        else:
            print("remove")
    
    clusterDF = pd.DataFrame(clusterList, columns=["gene", program])
    return(clusterDF)

######################################################################
# function 3: wrapper function for first 2 functions                 #
#             returns a dataframe with genes linked with the cluster #
######################################################################
def masterDFwrapper(filePath, OGprogram, columnOG):
    sortDF = sortGenes(filePath, OGprogram)
    if OGprogram == "orthofinder":
        sortDF = sortDF.drop(['OG', 'Gene Tree Parent Clade'], axis=1)

    cluster_df = clusterGeneList(sortDF, columnOG)
    return(cluster_df)

###################################################################################
# function 4: matching AT gene ID to the protein name, used in next few functions #
###################################################################################

def atFeature(atgene):
    at_features = "/path_to_ATfeatures/at_features.tsv"
    features = open(at_features)
    for line in features:
        genes = line.split("\t")
        if genes[0] == atgene.split(".")[0]:
            return genes[1] 

########################################################################################
# function 5: function to count number of clusters that match                          #
#             between the AT gene and the other species gene homolog                   #
#             for determine the composition of the conservative orthogroup composition #
########################################################################################

def countMatchingCluster(atRowInfo, spRowInfo):
    count = 0
    if atRowInfo[1] == spRowInfo[1]:
        count = count + 1
    if atRowInfo[2] == spRowInfo[2]:
        count = count + 1
    if atRowInfo[3] == spRowInfo[3]:
        count = count + 1
    if atRowInfo[4] == spRowInfo[4]:
        count = count + 1
    if atRowInfo[5] == spRowInfo[5]:
        count = count + 1
    if atRowInfo[6] == spRowInfo[6]:
        count = count + 1
   
    return(count)

#############################################################################################
# function 6: inclusive option - gene was included in the OG if found in one of the outputs #
#             for each gene subset - reorganize table to get species-specific composition   #
#             make output like that of original orthogroup table composition                #
#############################################################################################

def reformattedOG_atleast1(AT_ref_gene, OG_cluster):
    AT_new_OG = {}
    geneList = OG_cluster["gene"].tolist()
    refAT_feature = atFeature(AT_ref_gene)
    AT_new_OG[AT_ref_gene]={'ATfeature': refAT_feature,'refGene': AT_ref_gene,'Aar':[], 'Ath':[], 'Bra':[], 'Chi':[],'Cru':[],'Csa':[], 'Sal':[],'Tar':[]}
    for gene2 in geneList:
        if gene2.startswith("AT"):
            AT_new_OG[AT_ref_gene]['Ath'].append(gene2)
        elif gene2.startswith("CAR"):
            AT_new_OG[AT_ref_gene]['Chi'].append(gene2)
        elif gene2.startswith("Car"):
            AT_new_OG[AT_ref_gene]['Cru'].append(gene2)
        elif gene2.startswith("gene"):
            AT_new_OG[AT_ref_gene]['Tar'].append(gene2)
    return(AT_new_OG)

###########################################################################################
# function 7: stringent/conservative option - gene has to be found in 5/6 OG outputs      #
#             for each gene subset - reorganize table to get species-specific composition #
#             make output like that of original orthogroup table composition              #
###########################################################################################

def reformattedOG_common(AT_ref_gene, OG_cluster, atRow):
    AT_new_OG = {}
    refAT_feature = atFeature(AT_ref_gene)
    AT_new_OG[AT_ref_gene]={'ATfeature': refAT_feature,'refGene': AT_ref_gene,'Aar':[], 'Ath':[], 'Bra':[], 'Chi':[],'Cru':[],'Csa':[], 'Sal':[],'Tar':[]}
                
    for index, row in OG_cluster.iterrows():
        gene2 = row[0]
        if gene2.startswith("AT"):
            Ath_match = countMatchingCluster(atRow, row)
            if Ath_match >= 5:
                AT_new_OG[AT_ref_gene]['Ath'].append(gene2)
        elif gene2.startswith("CAR"):
            Chi_match = countMatchingCluster(atRow, row)
            if Chi_match >= 5:
                AT_new_OG[AT_ref_gene]['Chi'].append(gene2)
        elif gene2.startswith("Car"):
            Cru_match = countMatchingCluster(atRow, row)
            if Cru_match >= 5:
                AT_new_OG[AT_ref_gene]['Cru'].append(gene2)
        elif gene2.startswith("gene"):
            Tar_match = countMatchingCluster(atRow, row)
            if Tar_match >= 5:
                AT_new_OG[AT_ref_gene]['Tar'].append(gene2)
                
    return(AT_new_OG)
    
############################################################################
# function 8: to determine most inclusive OG composition for each scenario #
############################################################################
def maxSubsetGroup(subsetList, OG_df): 
    programMax = max(subsetList, key=lambda x: x[2])[0]
    clusterMax = max(subsetList, key=lambda x: x[2])[1]
    subset = OG_df.loc[OG_df[programMax]==clusterMax]
    subset = subset.fillna("*")
    return(subset)

################
################
## START HERE ##
################
################

#######################################################################
# step 1: read in files and call use wrapper function to get          #
#         from original OG files to dataframe, list, and dictionaries #
#######################################################################

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
BR_df = masterDFwrapper(broccoli, "broccoli", "BR")
OFb_df = masterDFwrapper(orthofinder_b, "orthofinder", "OFb")
OFd_df = masterDFwrapper(orthofinder_d, "orthofinder", "OFd")
OFm_df = masterDFwrapper(orthofinder_m, "orthofinder", "OFm")
ON_df= masterDFwrapper(orthnet, "orthnet", "ON")
SPd_df = masterDFwrapper(sonicparanoid_d, "sonicparanoid", "SPd")
SPm_df = masterDFwrapper(sonicparanoid_m, "sonicparanoid", "SPm")

# merging dataframes together
merge1 = BR_df.merge(OFb_df, on="gene", how="outer")
merge2 = merge1.merge(OFd_df, on="gene", how="outer")
merge3 = merge2.merge(OFm_df, on="gene", how="outer")
merge4 = merge3.merge(ON_df, on="gene", how="outer")
merge5 = merge4.merge(SPd_df, on="gene", how="outer")
merge6 = merge5.merge(SPm_df, on="gene", how="outer")

#######################################################
# step 2: calculate rand score, adjusted rand score,  # 
#         and jaccard index                           #
#######################################################
# calculate the rand score, adjusted rand score, and jaccard index
# include ON

metricsList=[]
merge6["SPd"] = merge6["SPd"].astype("str").replace('nan',np.nan)
merge6["SPm"] = merge6["SPm"].astype("str").replace('nan',np.nan)

for index, row in merge6.iterrows():
    gene = row[0]
    if gene.startswith("AT"):
        print(gene)
        atRow = row
    
        BR_cluster = row[1]
        OFb_cluster = row[2]
        OFd_cluster = row[3]
        OFm_cluster = row[4]
        ON_cluster = row[5]
        SPd_cluster = row[6]
        SPm_cluster = row[7]
        
        BR_subset = merge6.loc[merge6["BR"]==BR_cluster]
        OFb_subset = merge6.loc[merge6["OFb"]==OFb_cluster]
        OFd_subset = merge6.loc[merge6["OFd"]==OFd_cluster]
        OFm_subset = merge6.loc[merge6["OFm"]==OFm_cluster]
        ON_subset = merge6.loc[merge6["ON"]==ON_cluster]
        SPd_subset = merge6.loc[merge6["SPd"]==SPd_cluster]
        SPm_subset = merge6.loc[merge6["SPm"]==SPm_cluster]

        #calculate the rand score, adjusted rand score, proportion similarity
        programList = ["BR", "OFb", "OFd", "OFm", "SPd", "SPm", "ON"]
        clusterList = [BR_cluster, OFb_cluster, OFd_cluster, OFm_cluster, SPd_cluster, SPm_cluster, ON_cluster]
        subsetList = [BR_subset, OFb_subset, OFd_subset, OFm_subset, SPd_subset, SPm_subset, ON_subset]
        for i in range(len(programList)-1):
            for j in range(i+1,7):
                p1 = programList[i]
                p2 = programList[j]
                c1 = clusterList[i]
                c2 = clusterList[j]
                s1 = subsetList[i]
                s2 = subsetList[j]
                print(p1, p2)
                print(c1, c2)
                #print(s1, s2)

                # for getting correct sets for RS and ARS calculations
                # subset dataframe for each program; find all members in the same cluster; drop duplicates
                # duplicates - an issue for BR
                p1_df = merge6.loc[merge6[p1]==c1].drop_duplicates("gene")
                p2_df = merge6.loc[merge6[p2]==c2].drop_duplicates("gene")
                
                # join the two dataframes together
                union_df = pd.concat([p1_df,p2_df], ignore_index=True)
                union_df_nodup = union_df.drop_duplicates("gene")
                
                #print(union_df_nodup)
                union_df_nodup = union_df_nodup.fillna("*")

                # cluster values from each program's column - to use for calculating RS and ARS
                program1 = union_df_nodup[p1]
                program2 = union_df_nodup[p2]

                p1_clusters = list(set(union_df_nodup[p1].tolist()))
                p2_clusters = list(set(union_df_nodup[p2].tolist()))

                #print("program")
                #print(program1)
                #print(program2)
                
                # for calculating jaccard similarity index
                p1_only = set(merge6.loc[merge6[p1]==c1]["gene"].tolist())
                p2_only = set(merge6.loc[merge6[p2]==c2]["gene"].tolist())
                
                #print("p1, p2")
                #print(p1_only)
                #print(p2_only)
                
                # determine how many AT genes are in each orthogroup
                p1_AT = []
                for p1_gene in p1_only:
                    if p1_gene.startswith("AT"):
                        p1_AT.append(p1_gene)
                p1_AT_len = len(p1_AT)
                p1_AT_join = ",".join(p1_AT)
                
                p2_AT = []
                for p2_gene in p2_only:
                    if p2_gene.startswith("AT"):
                        p2_AT.append(p2_gene)
                p2_AT_len = len(p2_AT)
                p2_AT_join = ";".join(p2_AT)
            
                if len(p1_only) == 0 and len(p2_only) == 0:
                    jaccard = "na"
                    RS = "na"
                    ARS = "na"
                    intersect = "na"
                    un = "na"
                elif len(p1_only) == 0 or len(p2_only) == 0:
                    RS = "na"
                    ARS = "na"
                    intersect = len(p1_only & p2_only)
                    un = len(p1_only) + len(p2_only) - intersect
                    jaccard = float((intersect) / un)
                else:                
                    intersect = len(p1_only & p2_only)
                    un = len(p1_only) + len(p2_only) - intersect
                    jaccard = float((intersect) / un)
                    
                    RS = rand_score(program1, program2)
                    ARS = adjusted_rand_score(program1, program2)
               
                metricsList.append([gene,p1, c1, p2, c2, RS,ARS,intersect,un,jaccard,p1_AT_len,p1_AT_join,p2_AT_len,p2_AT_join])
                print(RS, ARS, intersect, un, jaccard, p1_AT_len, p2_AT_len)   

####################################### 
# step 3: convert list into dataframe #
#######################################
metricColumns = ["gene", "program1", "OG1", "program2", "OG2", "RS","ARS","intersection", "union", "Jaccard","p1_num_AT","p1_AT_genes","p2_num_AT","p2_AT_genes"]
metricsDF = pd.DataFrame(metricsList, columns = metricColumns)
          
metrics_out = main_path + "/output/diploid_pairwiseMetrics_240217.csv"
metricsDF.to_csv(metrics_out, sep = ",", index="False")  
    