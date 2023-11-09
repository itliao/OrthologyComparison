# 230703 - update paths 
import pandas as pd

# from: https://colab.research.google.com/github/zaneveld/full_spectrum_bioinformatics/blob/master/content/06_biological_sequences/reading_and_writing_fasta_files.ipynb#scrollTo=ZkpJZr86_9Jz
def fastaDictionary(FastaFile):
    fasta = open(FastaFile, "r")
    """Return a dict of {id:gene_seq} pairs based on the sequences in the input FASTA file
    input_file -- a file handle for an input fasta file
    """
    parsed_seqs = {}
    curr_seq_id = None
    curr_seq = []

    for line in fasta:
        line = line.strip()
        
        #set up the dictionary
        if line.startswith(">"):
            # before going on to the next sequence in the file, join the final sequences together in the dictionary
            if curr_seq_id is not None:
                parsed_seqs[curr_seq_id] = ''.join(curr_seq)
                
            # rename/update/ curr_seq_id if there is a new line with ">"
            curr_seq_id = line[1:]
            # reset the sequence list to join
            curr_seq = []
            continue
        
        #add the sequence to the list
        curr_seq.append(line)

    #Add the final sequence to the dict
    parsed_seqs[curr_seq_id] = ''.join(curr_seq)
    return(parsed_seqs)

def convertFastaChi(FastaFile, geneListFile, outputFasta):
    fastaDict = fastaDictionary(FastaFile)
    for gene_id, seq in fastaDict.items():
        #print(gene_id)
        #print(seq)
        primaryList = open(geneListFile, "r")
        for gene in primaryList:
            primaryTranscript = gene.strip("\n")
            #print(primaryTranscript)
            #print("gene: ", gene_id)
            if primaryTranscript == gene_id:
                #print(gene_id)
                newFasta = open(outputFasta, "a+")
                newFasta.write(">"+gene_id)
                newFasta.write("\n")
                newFasta.write(seq)
                newFasta.write("\n")
                newFasta.close()

def convertFastaCsa(FastaFile, geneListFile, outputFasta):
    fastaDict = fastaDictionary(FastaFile)
    for gene_id, seq in fastaDict.items():
        #print(gene_id)
        #print(seq)
        geneIDlist = gene_id.strip().split(" ")
        fullName = geneIDlist[0].split(".")
        transcriptName = fullName[0] + "."+ fullName[1]
        primaryList = open(geneListFile, "r")
        for gene in primaryList:
            primaryTranscript = gene.strip()
            #print(primaryTranscript)
            if primaryTranscript == transcriptName:
                newFasta = open(outputFasta, "a+")
                newFasta.write(">"+gene_id)
                newFasta.write("\n")
                newFasta.write(seq)
                newFasta.write("\n")
                newFasta.close()
                
main_path = r"/path_to_files"

#for Cardamine hirsuta - on the cluster
ChiFastaFile = main_path + "/Phytozome_Transcripts/carhr38.aa.fa"
ChiPrimaryTranscripts = main_path + "/Chi_primaryGeneID.txt"
ChiNewFasta_path = main_path + "/Chirsuta.primary.pep.fa"
convertFastaChi(ChiFastaFile, ChiPrimaryTranscripts, ChiNewFasta_path)

#for Camelina sativa - on the cluster
CsaFastaFile = main_path + "/Phytozome_Transcripts/Camelina_sativa.Cs.pep.all.fa"
CsaPrimaryTranscripts = main_path + "/Csa_primaryGeneID.txt"
CsaNewFasta_path = main_path + "/Csativa.primary.pep.fa"
convertFastaCsa(CsaFastaFile, CsaPrimaryTranscripts, CsaNewFasta_path)
