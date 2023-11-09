import pandas as pd
import numpy as np
import os
import re

ortho_path = r"/path_to_output/Tar2"
genome_path = r"/path_to_genome"

#path to files
protein_path = genome_path + "/Tarvense_v2_protein.fa" 
gff_path = genome_path + "/Tarvense_v2_genomic.gff" 
orthnet_gtf_path = ortho_path + "/Tar.gtfParsed.txt" 

# create a dataframe/file that has the gene id and the protein id
# input - original gff 

# 1) create list of lists of gene and protein id
geneProteinList = []

gff = open(gff_path, "r")
for line in gff:
    if not line.startswith("#"):
        newline = line.strip().split("\t")
        if newline[2] == "CDS":
            #print(newline[2])
            id_info = newline[-1]
            #print(id_info)
            list_info = id_info.split(";")
            for item in list_info:
                new_item = item.split("=")
                if new_item[0] == "Parent":
                    part_parent = new_item[1].split("-")
                    gene_id = "gene-"+part_parent[1]
                    #print(gene_id)
                elif new_item[0] == "protein_id":
                    protein_id = new_item[1]
                    #print(protein_id)
            
                    geneProteinList.append([gene_id, protein_id])
                
                
# 2) create dataframe and write to file
geneProteinDF = pd.DataFrame(geneProteinList, columns=["gene_id", "protein_id"])
geneProteinDF_noDup = geneProteinDF.drop_duplicates().sort_values(by="gene_id")
#geneProteinDF_noDup.to_csv(ortho_path+"/Tarvense_v2_geneProtein_pair.tsv", sep="\t", index=False)

# function for converting gene_id to protein_id 
def Gene2Protein(gene):
    for index, row in geneProteinDF_noDup.iterrows():
        if row[0] == gene:
            protein_id = row[1]
            return protein_id   

# function for converting protein_id to gene_id
def Protein2Gene(protein):
    for index, row in geneProteinDF_noDup.iterrows():
        if row[1] == protein:
            gene_id = row[0]
            return gene_id 
        
# convert protein fasta to have gene_id - takes a while, so may want to run this on the cluster.
fasta_out_path = ortho_path+"/Tar_v2_protein_geneID.fa"
fasta = open(protein_path,"r")
for line in fasta:
    if line.startswith(">"):
        newline = line.strip().split(" ")
        protein_id = newline[0].replace(">","")
        gene_id = Protein2Gene(protein_id)
        #print(protein_id)
        #print(gene_id)
        
        fasta_out = open(fasta_out_path, "a+")
        fasta_out.write(">"+gene_id+"\n")
        fasta_out.close()
    else:
        fasta_out = open(fasta_out_path, "a+")
        fasta_out.write(line)
        fasta_out.close()

# convert OrthNet gtf to have protein_id 
gtf_out_path = ortho_path+"/Tar_v2_proteinID.gtfParsed.txt"
gtf = open(orthnet_gtf_path,"r")
for line in gtf:
    if line.startswith("gene-"):
        newline = line.strip().split("\t")
        gene_id = newline[0]
        protein_id = Gene2Protein(gene_id)
        #print(protein_id)
        #print(gene_id)
        
        gtf_out = open(gtf_out_path, "a+")
        gtf_out.write(protein_id+"\t")

        for i in range(len(newline)):
            #print(i)
            if i < len(newline)-2:
                gtf_out.write(newline[i+1]+"\t")
            elif i == len(newline)-1:
                gtf_out.write(newline[i]+"\n")
                gtf_out.close()
    else:
        gtf_out = open(gtf_out_path, "a+")
        gtf_out.write(line)
        gtf_out.close()
