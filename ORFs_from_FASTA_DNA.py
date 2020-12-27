#!/usr/bin/env python

# Enter file name as 'filename.txt'
def ORFs(filenamein, filenameout):
        sss = open(filenamein).read()
        f = open(filenameout, "w")
        
        # Parse the FASTA file into a dictionary
        global FASTA_seq
        FASTA_seq = {}
        strings = sss.strip().split('>')
        for s in strings:
            if len(s) == 0:
                continue
            parts = s.split()
            label = parts[0]
            bases = ''.join(parts[1:])
            FASTA_seq[label] = bases
        
        global results_dict
        results_dict = {}


        # DNA codon dictionary
        codon_dict_DNA = {'TTT': 'F', 'CTT': 'L', 'ATT': 'I', 'GTT': 'V',
                          'TTC': 'F', 'CTC': 'L', 'ATC': 'I', 'GTC': 'V',
                          'TTA': 'L', 'CTA': 'L', 'ATA': 'I', 'GTA': 'V',
                          'TTG': 'L', 'CTG': 'L', 'ATG': 'M', 'GTG': 'V',
                          'TCT': 'S', 'CCT': 'P', 'ACT': 'T', 'GCT': 'A',
                          'TCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
                          'TCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
                          'TCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
                          'TAT': 'Y', 'CAT': 'H', 'AAT': 'N', 'GAT': 'D',
                          'TAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
                          'TAA': 'Stop', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
                          'TAG': 'Stop', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
                          'TGT': 'C', 'CGT': 'R', 'AGT': 'S', 'GGT': 'G',
                          'TGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
                          'TGA': 'Stop', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
                          'TGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}

        
        # For each sequence in the dictionary
        for key in FASTA_seq:      
            # Create a variable that will contain all of the ORFs for this
            # sequence 
            list_of_results = []
            # Create a variable for the forward sequence
            forward = FASTA_seq[key]
            # Creat a variable for the reverse sequence
            reverse = forward[::-1].translate(str.maketrans('ACGT', 'TGCA'))
 
            # For each position in the forward sequence
            for i in range(0, len(forward), 1):
                # If there are at least three nucleotides left in the sequence
                if (len(forward) - i) >= 3:
                    # If the position is part of a "Start" codon
                    if forward[i:i+3] == 'ATG':
                        # Create the variable that will contain the ORF sequence
                        result = ''
                        # For this part of the sequence through the end:
                        for w in range(i, len(forward), 3):
                            # If there are still at least three nucleotides left
                            # (Not sure why this isn't redundant, but you need it)
                            if (len(forward) - w) >= 3:
                                # Find the amino acid symbol that matches
                                # the DNA codon
                                symbol = codon_dict_DNA[forward[w:w+3]]
                                # If it's not a stop codon, add the amino acid
                                # to protein sequence
                                if symbol != 'Stop':
                                    result += symbol
                                # If it is a stop codon, add the whole ORF
                                # sequence to the list of results
                                else:
                                    list_of_results.append(result)
                                    break
                                            
                # Same deal but for the reverse sequence
                if (len(reverse) - i) >= 3:    
                    if reverse[i:i+3] == 'ATG':
                        rev_result = ''
                        for z in range(i, len(reverse), 3):
                            if (len(reverse) - z) >= 3:
                                rev_symbol = codon_dict_DNA[reverse[z:z+3]]
                                if rev_symbol != 'Stop':
                                    rev_result += rev_symbol
                                else:
                                    list_of_results.append(rev_result)
                                    break

                # Add all of the individual ORFs for a single sequence to
                # a list
                ORFs_list = []
                for i in range(0, len(list_of_results), 1):
                    orf = list_of_results[i]
                    ORFs_list.append(orf)    
                # Remove duplicates from that list
                nonred_ORFs = list(set(ORFs_list))
                # Attach the sequence name to the list of ORFs in a dictionary
                results_dict[key] = nonred_ORFs
                        
            # At the very end, print the dictionary
            print (results_dict)
            
            # Write a file containing the sequence name followed by each
            # ORF on a new line
            f.write("> %s\n" %key)
            for i in results_dict[key]:
                f.write("%s\n" %i)
            f.write("\n")
        f.close()    
