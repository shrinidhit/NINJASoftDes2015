# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 11:24:42 2014

@author: YOUR NAME HERE

"""

# you may find it useful to import these variables (although you are not required to use them)
from amino_acids import aa, codons, aa_table, matchings
import random
from load import load_seq

def shuffle_string(s):
    """ Shuffles the characters in the input string
        NOTE: this is a helper function, you do not have to modify this in any way """
    return ''.join(random.sample(s,len(s)))

### YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    return matchings[nucleotide]

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
    
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # TODO: implement this
    complement =  "".join(get_complement(nucl) for nucl in dna)
    reverse_complement = complement[::-1]
    return reverse_complement

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.
        
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    ORF = [] #initializing ORF, a list of codons in reading frame
    for i in range(0, len(dna) -2, 3): #start indexes of all codons in DNA
        codon = dna[i:i+3] #current codon
        if aa_table[codon] != "|":
            ORF.append(codon)
        else: #if we have hit a stop codon
            return ''.join(ORF) #converts to string
    return dna #if no stop codon


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence and returns
        them as a list.  This function should only find ORFs that are in the default
        frame of the sequence (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG']
    """
    #loop through DNA with a checker i that goes by 3s. 
    #If we hit a start codon, append call rest_of_ORFs to list of frames
    #increase checker by length of appended ORF
    allORFS = []
    i = 0
    while i <= len(dna) - 3:
        codon = dna[i:i+3]
        if codon == 'ATG':
            ORF = rest_of_ORF(dna[i:])
            allORFS.append(ORF)
            i += len(ORF)
        else:
            i += 3
    return allORFS


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in the same frame,
        it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    allORFS = []
    allORFS.extend(find_all_ORFs_oneframe(dna))
    allORFS.extend(find_all_ORFs_oneframe(dna[1:]))
    allORFS.extend(find_all_ORFs_oneframe(dna[2:]))
    return allORFS

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    strand1 = dna
    strand2 = get_reverse_complement(dna)

    allORFS = []
    allORFS.extend(find_all_ORFs(strand1))
    allORFS.extend(find_all_ORFs(strand2))

    return allORFS


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
        Also returns multiple ORFS if multiple have the same length
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    possibleORFS = find_all_ORFs_both_strands(dna)
    longORF = possibleORFS[0]
    length = len(longORF)
    for ORF in possibleORFS:
        if length < len(ORF):
            longORF = ORF
            length = len(ORF)
        elif length == len(ORF):
            longORF = [longORF, ORF]
    return longORF

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF 
    """

    ORFLengths = []
    for i in range(0,num_trials):
        testdna = shuffle_string(dna)
        if type(longest_ORF) == str:
            ORFLengths.append(len(longest_ORF))
        else:
            ORFLengths.append(len(longest_ORF[0]))
    return max(ORFLengths)

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).
        
        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    proteins = []
    for i in range(0, len(dna)-2, 3):
        aa = dna[i:i+3]
        proteins.append(aa_table[aa])
    return ''.join(proteins)


def gene_finder(dna, threshold):
    """ Returns the amino acid sequences coded by all genes that have an ORF
        larger than the specified threshold.
        
        dna: a DNA sequence
        threshold: the minimum length of the ORF for it to be considered a valid
                   gene.
        returns: a list of all amino acid sequences whose ORFs meet the minimum
                 length specified.
    """
    # TODO: implement this
    longORFs = []
    ORFs = find_all_ORFs_both_strands(dna)
    for ORF in ORFs:
        if len(ORF) >= threshold:
            longORFs.append(coding_strand_to_AA(ORF))
    return longORFs


if __name__ == "__main__":
    import doctest
    doctest.testmod()