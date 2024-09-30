def transcription(dna_seq):
    complement_dic = {'T': 'A', 'A': 'T', 'G': 'C', 'C': 'G'}
    complement_dna = ""

    # create the complement dna and mrna from the dna seq
    for nucleotide in dna_seq.upper():
        if nucleotide in complement_dic:
            complement_dna += complement_dic.get(nucleotide)
            
    mrna = complement_dna.replace('T', 'U')

    return complement_dna, mrna

def translate_codon(mrna):
    mrna_codon_dict = {
        'UUU': ['Phe', 'F'], 'UUC': ['Phe', 'F'],               # Phenylalanine
        'UUA': ['Leu', 'L'], 'UUG': ['Leu', 'L'],               # Leucine
        'CUU': ['Leu', 'L'], 'CUC': ['Leu', 'L'], 
        'CUA': ['Leu', 'L'], 'CUG': ['Leu', 'L'],               # Leucine (more codons)
        'AUU': ['Ile', 'I'], 'AUC': ['Ile', 'I'], 'AUA': ['Ile', 'I'],  # Isoleucine
        'AUG': ['Met', 'M'],                                    # Methionine (start codon)
        'GUU': ['Val', 'V'], 'GUC': ['Val', 'V'], 'GUA': ['Val', 'V'], 'GUG': ['Val', 'V'],  # Valine
        'UCU': ['Ser', 'S'], 'UCC': ['Ser', 'S'], 'UCA': ['Ser', 'S'], 'UCG': ['Ser', 'S'],  # Serine
        'AGU': ['Ser', 'S'], 'AGC': ['Ser', 'S'],               # Serine (more codons)
        'CCU': ['Pro', 'P'], 'CCC': ['Pro', 'P'], 'CCA': ['Pro', 'P'], 'CCG': ['Pro', 'P'],  # Proline
        'ACU': ['Thr', 'T'], 'ACC': ['Thr', 'T'], 'ACA': ['Thr', 'T'], 'ACG': ['Thr', 'T'],  # Threonine
        'GCU': ['Ala', 'A'], 'GCC': ['Ala', 'A'], 'GCA': ['Ala', 'A'], 'GCG': ['Ala', 'A'],  # Alanine
        'UAU': ['Tyr', 'Y'], 'UAC': ['Tyr', 'Y'],               # Tyrosine
        'UAA': ['Stop', 'Stop'], 'UAG': ['Stop', 'Stop'], 'UGA': ['Stop', 'Stop'],  # Stop codons
        'CAU': ['His', 'H'], 'CAC': ['His', 'H'],               # Histidine
        'CAA': ['Gln', 'Q'], 'CAG': ['Gln', 'Q'],               # Glutamine
        'AAU': ['Asn', 'N'], 'AAC': ['Asn', 'N'],               # Asparagine
        'AAA': ['Lys', 'K'], 'AAG': ['Lys', 'K'],               # Lysine
        'GAU': ['Asp', 'D'], 'GAC': ['Asp', 'D'],               # Aspartic Acid
        'GAA': ['Glu', 'E'], 'GAG': ['Glu', 'E'],               # Glutamic Acid
        'UGU': ['Cys', 'C'], 'UGC': ['Cys', 'C'],               # Cysteine
        'UGG': ['Trp', 'W'],                                    # Tryptophan
        'CGU': ['Arg', 'R'], 'CGC': ['Arg', 'R'], 'CGA': ['Arg', 'R'], 'CGG': ['Arg', 'R'],  # Arginine
        'AGA': ['Arg', 'R'], 'AGG': ['Arg', 'R'],               # Arginine (more codons)
        'GGU': ['Gly', 'G'], 'GGC': ['Gly', 'G'], 'GGA': ['Gly', 'G'], 'GGG': ['Gly', 'G']   # Glycine
    }

    amino_acids = []

    # get the codon by splitting the mrna into 3 for each codon
    for i in range(0, len(mrna), 3):
        codon = mrna[i:i+3]  # extract a 3-letter codon
        amino_acid_info = mrna_codon_dict.get(codon)

        if amino_acid_info:
            amino_acid, symbol = amino_acid_info
            amino_acids.append(f'{amino_acid} ({symbol})')
        else:
            amino_acids.append(f'Invalid codon: {codon}')

    return ' - '.join(amino_acids)

dna_seq = input("Input DNA = ")

while len(dna_seq) % 3 != 0:
    print("\nDNA should be multiple of 3")
    dna_seq = input("Reinput DNA = ")

dna_complement, mrna = transcription(dna_seq)
amino_acid_sequence = translate_codon(mrna)

print(f'\nComplement = {dna_complement}')
print(f'mRNA = {mrna}')
print(f'Aminoacid = {amino_acid_sequence}')