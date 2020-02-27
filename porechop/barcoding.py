from collections import defaultdict, Counter
from porechop.misc import load_fasta, print_table
from porechop.adapters import Adapter
import random


def makeKmers(fasta_file, k):
    """ 
    This function builds default dict of all k-mers in a fasta file
    """
    fastaSeqs = load_fasta(fasta_file)
    singleSeqDict = defaultdict(list)
    
    for seq in fastaSeqs:
        for base in range(len(seq[1])):
            #Generates Unique barcodes
            if len((seq[1][base:base+k])) >= k and seq[1][base:base+k] not in singleSeqDict.values():
                singleSeqDict[seq[0]].append(seq[1][base:base+k].upper())           
        
    return singleSeqDict
   
def reverseComplement(seq):
    """
    Get the reverse complement of a nucleotide sequence
    """
    
    # Translate bases to comp
    compStrDNA = str.maketrans('ACGTacgt', 'TGCAtcga')

    # Translate then get the reverse
    return seq.translate(compStrDNA)[::-1]

def findUnique(fastaFile, k):
    '''
    Takes a fasta file containing plasmid sequences to find unique sequences to use as barcodes for Porechop.
    Kmers of length K are generated for each sequences and then using set interection it returns a defaultDict
    of type set of each plasmid and a set of kmers that are unique to each plasmid.
    '''
    
    kmerList = makeKmers(fastaFile, k)  
    uniqueSetDict  = defaultdict(set)
    for plasmid in kmerList.keys():
        #For each plasmid
        differenceList = []
        #Build a list of barcodes from all other plasmids
        for i in kmerList:
            if i != plasmid:
                [differenceList.append(k) for k in kmerList[i]]
        
        #Enters a set made of barcodes from this plasmid only into a dict at the key of the plasmid
        uniqueSetDict[plasmid] = set(kmerList[plasmid]).difference(set(differenceList))
    
    return uniqueSetDict


# def filterBarcodes(uniqueSetDict):
#     '''
#     Similar to findUnique, this function goes through the unique set for each plasmid
#     and filters based on hamming distance. 
#     '''
#     filterBarcodes = defaultdict(list)
    
#     for plasmid in uniqueSetDict:
#         compareList = []
#         filterBarcodes[plasmid] = []
#         for otherPlasmid in uniqueSetDict:
#             if plasmid != otherPlasmid:
#                 [compareList.append(i) for i in uniqueSetDict[otherPlasmid]]
                
#         for barcode in uniqueSetDict[plasmid]:
#             [filterBarcodes[plasmid].append(barcode) for i in compareList if 4 < hamming_distance(barcode, i) < 1]
    
    
#     return filterBarcodes

def buildBarcodes(uniqueSetDict):
    
    pairing = defaultdict(tuple)

    for x in uniqueSetDict:
        
        forwardBC = random.choice(list(uniqueSetDict[x]))
        reverseBC = reverseComplement(forwardBC)
        pairing[x] = (forwardBC, reverseBC)
        
        
    return pairing

def formatBarcodes(fastaFile, k = 24):
    
    pairing = buildBarcodes(findUnique(fastaFile, k))
    
    formattedBarcodeList = [Adapter('Rapid', start_sequence=('Rapid_adapter',
                                                             'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA'))]
    N = 0
    
    for barcode in pairing:
        shortName  = 'BC'+str(N)
        formattedBarcodeList.append(Adapter('Barcode '+barcode+'(forward)',
                                            start_sequence = (shortName+'_rev' , pairing[barcode][0]), 
                                           end_sequence = (shortName, pairing[barcode][1])))
        
        formattedBarcodeList.append(Adapter('Barcode '+barcode+'(reverse)',
                                            start_sequence = (shortName, pairing[barcode][1]), 
                                           end_sequence = (shortName+'_rev', pairing[barcode][0])))
    
        N += 1
        
    
    return formattedBarcodeList

def hamming_distance(sequence1, sequence2): 
    ''' Function to calculate Hamming distance between two kmers.
    Hamming distance is defined as the total number of changes required
     to make the sequences identical, or alternatively the total number
     of bases that do not match between the sequences.
    '''
    distance = 0
    
    #Compares each position in both seqs and counts each position that is different
    for x in range(len(sequence1)):
        if sequence1[x] != sequence2[x]:
            distance += 1
        
    return distance


    