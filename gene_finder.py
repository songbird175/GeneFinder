# -*- coding: utf-8 -*-
"""
This package contains a set of functions which can be used to find genes within a given genome.
Results can be compared to those from the BLAST website.

Celina Bekins

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """

    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'


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
    return_list = [] #this list is what will eventually be turned into the final string

    for letter in dna:
        new_letter = get_complement(letter)
        return_list.append(new_letter) #adds each new compliment to the end of the list

    before = ''.join(return_list) #turns the list into a string again, before flipping it around
    after = before[::-1] #flips the string around (so it's backwards)

    return after


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    for index in range(0, len(dna), 3): #loops through, stopping every three indices
        if dna[index:index+3] == 'TAA' or dna[index:index+3] == 'TAG' or dna[index:index+3] == 'TGA':
            return dna[0:index] #returns up to, but not including, the stop codon
    return dna #if there's no stop codon, returns the string


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("CCCATGATGGACATATGATT")
    ['ATGATGGACATA']
    >>> find_all_ORFs_oneframe("ATGCATGGATGATT")
    ['ATGCATGGA']
    """
    ORF_list = []
    index = 0
    while index < len(dna): #because a for loop doesn't allow you to modify the index...
        if dna[index:index+3] == 'ATG':
            list_add = rest_of_ORF(dna[index:])
            ORF_list.append(list_add)
            index += len(list_add) #the string added to the list will have a length divisible by 3, so testing after it automatically stays w/in the frame
        else:
            index += 3 #this also makes sure that we stay within the original reference frame
    return ORF_list



def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    ref1 = find_all_ORFs_oneframe(dna) #start in the first reference frame
    ref2 = find_all_ORFs_oneframe(dna[1:]) #moves to the next one
    ref3 = find_all_ORFs_oneframe(dna[2:]) #and the third
    entire_list = ref1 + ref2 + ref3 #not a list of lists; instead, just lists added together to make one
    return entire_list


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    strand1 = find_all_ORFs(dna)
    strand2 = find_all_ORFs(get_reverse_complement(dna))
    both_strands = strand1 + strand2
    return both_strands


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    allORFs = find_all_ORFs_both_strands(dna)
    first = allORFs[0] #start with assuming that the first item in the list is the largest
    for item in allORFs:
        if len(first) <= len(item):
            first = item #updates the largest term to be the next one if it's larger than the current largest
    return first #yay


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    import random
    first = 'ATG' #initializing the first string as a random string
    for num in range(num_trials): #test until we've reached num_trials tests
        test = ''.join(random.sample(dna, len(dna))) #this shuffles the string and puts it back into a single string again
        longer = longest_ORF(test) #here we add the longest ORF found in the shuffled sequence to a list of long ORFs
        if len(first) <= len(longer): #comparing lengths just like in longest_ORF
            first = longer
    return len(first)


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents a protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    import amino_acids
    proteinList = []
    for index in range(0, len(dna), 3): #runs through in multiples of three to get full codons
        if len(dna) - index >= 3:
            AminoAcid = amino_acids.aa_table[dna[index:index+3]] #gets the amino acid associated with each codon
            proteinList.append(AminoAcid)
    final = ''.join(proteinList) #puts the list together as a single string
    return final


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    list_aminos = []
    threshold = longest_ORF_noncoding(dna, 1500)
    actualORFs = find_all_ORFs_both_strands(dna) #these are all the reference frames we want to look at
    for item in actualORFs:
        if len(item) > threshold: #if the reference frame is larger than the threshold length we established above using longest_ORF_noncoding...
            list_aminos.append(coding_strand_to_AA(item)) #we compute the amino acids for that reference frame and add that to a list
    return list_aminos #FINALLY

if __name__ == "__main__":
    import doctest
    doctest.testmod()

from load import load_seq
dna = load_seq("./data/X73525.fa")

print gene_finder(dna)
