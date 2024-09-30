import itertools

def generate_mrna_combinations_and_count(aminoacid):
    encode_aminoacid_dict = {
        'F': ['UUU', 'UUC'],                                       # Phenylalanine
        'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],           # Leucine
        'I': ['AUU', 'AUC', 'AUA'],                                # Isoleucine
        'M': ['AUG'],                                              # Methionine (start codon)
        'V': ['GUU', 'GUC', 'GUA', 'GUG'],                         # Valine
        'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],           # Serine
        'P': ['CCU', 'CCC', 'CCA', 'CCG'],                         # Proline
        'T': ['ACU', 'ACC', 'ACA', 'ACG'],                         # Threonine
        'A': ['GCU', 'GCC', 'GCA', 'GCG'],                         # Alanine
        'Y': ['UAU', 'UAC'],                                       # Tyrosine
        'Stop': ['UAA', 'UAG', 'UGA'],                             # Stop codons
        'H': ['CAU', 'CAC'],                                       # Histidine
        'Q': ['CAA', 'CAG'],                                       # Glutamine
        'N': ['AAU', 'AAC'],                                       # Asparagine
        'K': ['AAA', 'AAG'],                                       # Lysine
        'D': ['GAU', 'GAC'],                                       # Aspartic Acid
        'E': ['GAA', 'GAG'],                                       # Glutamic Acid
        'C': ['UGU', 'UGC'],                                       # Cysteine
        'W': ['UGG'],                                              # Tryptophan
        'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],           # Arginine
        'G': ['GGU', 'GGC', 'GGA', 'GGG']                          # Glycine
    }

    # generate a list of possible codon combinations for each amino acid
    codon_combinations = [encode_aminoacid_dict[aa] for aa in aminoacid]  

    # generate all possible combinations of codons
    all_combinations = list(itertools.product(*codon_combinations))

    # for each mrna seq, calculate the frequencies of the codon
    for combination in all_combinations:
        mrna_sequence = "".join(combination)
        codon_count = {}

        for codon in combination:
            if codon in codon_count:
                codon_count[codon] += 1
            else:
                codon_count[codon] = 1

        print("mRNA =", mrna_sequence)
        for codon, count in codon_count.items():
            print(f"{codon} = {count}")
        print("\n" + "-"*30 + "\n")

aminoacid = input("Input Aminoacid = ")

while len(aminoacid) > 3:
    print("\nMax 3 aminoacids input ")
    aminoacid = input("Reinput Aminoacid = ")

generate_mrna_combinations_and_count(aminoacid)