NA_BASES = 'ATCGU'
DNA_BASES = 'ATCG'
RNA_BASES = 'AUCG'
stop_code = '_'
codon_lenght = 3

genetic_code = {
    'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': '*', 'ACA': 'U', 'ACC': 'T', 'ACG': 'U', 'ACU': 'U',
    'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
    'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
    'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
    'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
    'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
    'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S', 'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
    'UAC': 'Y', 'UAU': 'Y', 'UAA': '_', 'UAG': '_', 'UGC': 'C', 'UGU': 'C', 'UGA': '_', 'UGG': 'W',
}

tableaaweights = {  # weight of residues, i.e. weight of water already removed
    'A': 071.09,  # alanine
    'R': 156.19,  # arginine
    'D': 114.11,  # aspartic acid
    'N': 115.09,  # asparagine
    'C': 103.15,  # cysteine
    'E': 129.12,  # glutamic acid
    'Q': 128.14,  # glutamine
    'G': 057.05,  # glycine
    'H': 137.14,  # histidine
    'I': 113.16,  # isoleucine
    'L': 113.16,  # leucine
    'K': 128.17,  # lysine
    '*': 131.19,  # methionine
    'M': 131.19,  # methionine
    'F': 147.18,  # phenylalanine
    'P': 097.12,  # proline
    'S': 087.08,  # serine
    'T': 101.11,  # threonine
    'W': 186.12,  # tryptophan
    'Y': 163.18,  # tyrosine
    'V': 099.14,  # valine
}

DNA_nucleotides_complements = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G',
}

RNA_nucleotides_complements = {
    'A': 'U',
    'U': 'A',
    'G': 'C',
    'C': 'G',
}


def is_only_characters(sequence: str, allowed_characters: str) -> bool:
    for character in sequence:
        if character not in allowed_characters:
            return False
    return True


def is_valid_na(sequence: str) -> bool:
    return is_only_characters(sequence, allowed_characters=NA_BASES)


def is_valid_dna(sequence: str) -> bool:
    return is_only_characters(sequence, allowed_characters=DNA_BASES)


def is_valid_rna(sequence: str) -> bool:
    return is_only_characters(sequence, allowed_characters=RNA_BASES)


def keep_only_characters(sequence: str, allowed_characters: str) -> str:
    for character in sequence:
        if character not in allowed_characters:
            sequence = sequence.replace(character, '')
    return sequence


def reverse_sequence(sequence: str) -> str:
    return sequence[::-1]


def transform_dnabases_to_rna(dna_sequence: str) -> str:
    return dna_sequence.replace('T', 'U')


def transform_rnabases_to_dna(rna_sequence: str) -> str:
    return rna_sequence.replace('U', 'T')


def complement_sequence(sequence: str, rna_out: bool = True) -> str:
    complemented_sequence = ""
    if rna_out is True:
        converted_in_rna_bases_sequence = transform_dnabases_to_rna(sequence)
        for nucleotide in converted_in_rna_bases_sequence:
            complemented_sequence += RNA_nucleotides_complements.get(nucleotide, '?')
    else:
        converted_in_dna_bases_sequence = transform_rnabases_to_dna(sequence)
        for nucleotide in converted_in_dna_bases_sequence:
            complemented_sequence += DNA_nucleotides_complements.get(nucleotide, '?')
    return complemented_sequence


def correct_and_complement_sequence(sequence: str) -> str:
    corrected_and_complemented_sequence = ""
    if is_valid_na is False:
        keep_only_characters(sequence, allowed_characters=NA_BASES)
        corrected_and_complemented_sequence = complement_sequence(corrected_sequence)
    return corrected_and_complemented_sequence


def get_nucleotide_count(sequence: str) -> {}:
    nucleotide_count = {'A': sequence.count('A'), 'T': sequence.count('T'), 'C': sequence.count('C'),
                        'G': sequence.count('G')}
    return nucleotide_count


def build_dna_stats(sequence: str) -> {}:
    sequence_length = len(sequence)
    nucleotide_count = get_nucleotide_count(sequence)
    nucleotide_relative_abundance = {'A': nucleotide_count.get('A') / sequence_length,
                                     'T': nucleotide_count.get('T') / sequence_length,
                                     'C': nucleotide_count.get('C') / sequence_length,
                                     'G': nucleotide_count.get('G') / sequence_length}
    return nucleotide_relative_abundance


def compute_oligonucleotide_melt_temperature(sequence: str) -> float:
    nucleotide_count = get_nucleotide_count(sequence)
    if len(sequence) < 14:
        oligonucleotide_melt_temperature = (nucleotide_count.get('A') + nucleotide_count.get('T')) * 2 + (
                nucleotide_count.get('C') + nucleotide_count.get('G')) * 4
    else:
        oligonucleotide_melt_temperature = 64.9 + 41 * \
                                           (nucleotide_count.get('G') + nucleotide_count.get('C') - 16.4) \
                                           / (nucleotide_count.get('A') + nucleotide_count.get('T')
                                              + nucleotide_count.get('G')
                                              + nucleotide_count.get('C'))
    return round(oligonucleotide_melt_temperature, 1)


def split_sequence_to_groups_of_characters(sequence: str, characters_group_length: int,
                                           frame_starting_character_position=1) -> []:
    characters_groups = []
    for character_index in range(
            frame_starting_character_position - 1, len(sequence), characters_group_length):
        characters_groups.append(sequence[character_index: character_index + characters_group_length])
    return characters_groups


def split_na_sequence_to_codons(sequence: str, frame: int = 1) -> []:
    codons = split_sequence_to_groups_of_characters(
        sequence, characters_group_length=codon_lenght, frame_starting_character_position=frame
    )
    if len(codons[-1]) != codon_lenght:
        del codons[-1]  # le dernier codon doit avoir une longueur de 3 pour être dans la liste des codons
    return codons


def prepare_column_print(sequence: str, maximum_number_of_characters_per_column: int = 10,
                         maximum_number_of_columns_per_line: int = 6) -> str:
    columns_to_print = ''
    columns_of_characters_in_the_line = []
    lines = []
    maximum_number_of_characters_per_line = \
        maximum_number_of_characters_per_column * maximum_number_of_columns_per_line
    characters_groups = split_sequence_to_groups_of_characters(
        sequence, maximum_number_of_characters_per_column)
    first_character_of_the_line_number = 1
    for character_group in characters_groups:
        columns_of_characters_in_the_line.append(character_group)
        if len(columns_of_characters_in_the_line) == maximum_number_of_columns_per_line:
            lines.append(columns_of_characters_in_the_line)
            columns_of_characters_in_the_line = []
    if len(
            columns_of_characters_in_the_line) < maximum_number_of_characters_per_line \
            and columns_of_characters_in_the_line != []:
        lines.append(columns_of_characters_in_the_line)
    for line in lines:
        columns_to_print += '{:0{padding_width}} '.format(
            first_character_of_the_line_number, padding_width=len(str(len(sequence)))
        )
        columns_to_print += ' '.join(line)
        if line != lines[-1]:
            columns_to_print += '\n'
        first_character_of_the_line_number += maximum_number_of_characters_per_line
    return columns_to_print


def translate_rna_to_peptide(rna_sequence: str, stopatstopcode: bool = True):
    peptide = ""
    codons = split_na_sequence_to_codons(rna_sequence)
    if stopatstopcode is False:
        for codon in codons:
            translated_codon = genetic_code.get(codon, '?')
            peptide += translated_codon
    else:
        for codon in codons:
            translated_codon = genetic_code.get(codon, '?')
            if translated_codon == stop_code:
                break
            peptide += translated_codon

    return peptide


def get_peptide_weight(peptide: str) -> float:
    water_weight = 18.00
    peptide_weight = 0.00
    for aa in peptide:
        peptide_weight += tableaaweights.get(aa, 0)
    peptide_weight += water_weight
    return round(peptide_weight, 2)


nucleotides_sequence = str(input("Entrez une séquence de nucléotides:"))
if is_valid_na(nucleotides_sequence) is False:
    corrected_sequence = keep_only_characters(nucleotides_sequence, allowed_characters=NA_BASES)
if is_valid_dna(nucleotides_sequence) is True:
    dna = nucleotides_sequence
    rna = complement_sequence(reverse_sequence(dna), rna_out=True)
else:
    rna = nucleotides_sequence
    dna = complement_sequence(reverse_sequence(rna), rna_out=False)
peptide = translate_rna_to_peptide(rna, stopatstopcode=False)
peptide_weight = get_peptide_weight(peptide)
dna_melt_temperature = compute_oligonucleotide_melt_temperature(dna)
dna_bases_relative_abundance = build_dna_stats(dna)
prepared_columns_dna = prepare_column_print(dna, maximum_number_of_characters_per_column=10,
                                            maximum_number_of_columns_per_line=6)
prepared_columns_rna = prepare_column_print(rna, maximum_number_of_characters_per_column=10,
                                            maximum_number_of_columns_per_line=6)
prepared_columns_peptide = prepare_column_print(peptide, maximum_number_of_characters_per_column=10,
                                                maximum_number_of_columns_per_line=6)
print()

print("DNA sequence: {}{}".format('\n', prepared_columns_dna))
print()
print("DNA stats: {}".format(dna_bases_relative_abundance))
print()
print("DNA meling point: {}".format(dna_melt_temperature))
print()
print("RNA sequence: {}{}".format('\n', prepared_columns_rna))
print()
print("Peptide: {}{}".format('\n', prepared_columns_peptide))
print()
print("Peptide weight: {}".format(peptide_weight))
