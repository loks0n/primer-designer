from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

def translate_dna_to_protein(dna_sequence):
    """Translate DNA sequence to protein sequence."""
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }

    protein = ""
    for i in range(0, len(dna_sequence), 3):
        if i + 3 <= len(dna_sequence):
            codon = dna_sequence[i:i+3].upper()
            if codon in codon_table:
                protein += codon_table[codon]
            else:
                protein += 'X'  # Unknown amino acid
    return protein

def reverse_complement(dna_sequence):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement.get(base, base) for base in reversed(dna_sequence))

def find_codon_position(protein_position, dna_sequence):
    """Find the position of the codon in the DNA sequence for a given protein position."""
    return (protein_position - 1) * 3

def calculate_gc_content(sequence):
    """Calculate the GC content of a DNA sequence as a percentage."""
    gc_count = sum(1 for base in sequence.upper() if base in 'GC')
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0

def get_codons_for_amino_acid(amino_acid):
    """Get all possible codons for a given amino acid."""
    codon_map = {
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'N': ['AAT', 'AAC'],
        'D': ['GAT', 'GAC'],
        'C': ['TGT', 'TGC'],
        'Q': ['CAA', 'CAG'],
        'E': ['GAA', 'GAG'],
        'G': ['GGT', 'GGC', 'GGA', 'GGG'],
        'H': ['CAT', 'CAC'],
        'I': ['ATT', 'ATC', 'ATA'],
        'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
        'K': ['AAA', 'AAG'],
        'M': ['ATG'],
        'F': ['TTT', 'TTC'],
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],
        'W': ['TGG'],
        'Y': ['TAT', 'TAC'],
        'V': ['GTT', 'GTC', 'GTA', 'GTG'],
        '*': ['TAA', 'TAG', 'TGA']
    }
    return codon_map.get(amino_acid.upper(), [])

def design_mutagenesis_primers(dna_sequence, mutations):
    """
    Design primers for site-directed mutagenesis.

    Parameters:
    dna_sequence (str): The GPCR DNA sequence
    mutations (list): List of tuples (original_aa, position, desired_aa)

    Returns:
    list: List of dictionaries containing mutation details and primers

    The primers are designed following these principles:
    1. The mutation site (codon change) is positioned in the middle of the primer
    2. Typically 15-18 nucleotides flanking each side of the mutation
    3. Only the target codon is changed, surrounding sequence remains identical
    """
    dna_sequence = dna_sequence.upper()
    protein_sequence = translate_dna_to_protein(dna_sequence)
    results = []

    for original_aa, position, desired_aa in mutations:
        # Verify the original amino acid
        if position <= 0 or position > len(protein_sequence) or (protein_sequence[position-1] != original_aa):
            # Create error message safely without referencing potentially out-of-range positions
            error_msg = "Original amino acid mismatch or position out of range."
            if position > 0 and position <= len(protein_sequence):
                error_msg = f"Original amino acid mismatch or position out of range. Expected {protein_sequence[position-1]} at position {position}."
            
            results.append({
                'mutation': f"{original_aa}{position}{desired_aa}",
                'error': error_msg
            })
            continue

        # Find the codon position in DNA
        codon_pos = find_codon_position(position, dna_sequence)
        original_codon = dna_sequence[codon_pos:codon_pos+3]

        # Generate all possible codons for the desired amino acid
        desired_codons = get_codons_for_amino_acid(desired_aa)

        best_primers = None
        best_tm_diff = float('inf')
        best_gc_balance = float('inf')

        for desired_codon in desired_codons:
            # Try different flank lengths from 15-18 bases
            for flank_length in range(15, 19):
                # Ensure we don't go out of bounds
                left_flank = dna_sequence[max(0, codon_pos - flank_length):codon_pos]
                right_flank = dna_sequence[codon_pos+3:min(len(dna_sequence), codon_pos + 3 + flank_length)]

                # Create primers with the mutation
                forward_primer = left_flank + desired_codon + right_flank
                reverse_primer = reverse_complement(forward_primer)

                # Calculate properties
                fw_tm = mt.Tm_NN(Seq(forward_primer))
                rv_tm = mt.Tm_NN(Seq(reverse_primer))
                fw_gc = calculate_gc_content(forward_primer)
                rv_gc = calculate_gc_content(reverse_primer)

                # Check if primers meet relaxed requirements
                tm_requirement = fw_tm >= 65 and rv_tm >= 65 and abs(fw_tm - rv_tm) <= 5
                gc_requirement = 35 <= fw_gc <= 70 and 35 <= rv_gc <= 70

                if tm_requirement and gc_requirement:
                    # Calculate how close these primers are to ideal
                    tm_diff = abs(fw_tm - 70) + abs(rv_tm - 70)
                    gc_balance = abs(50 - fw_gc) + abs(50 - rv_gc)

                    # Choose the best primers
                    if tm_diff < best_tm_diff or (tm_diff == best_tm_diff and gc_balance < best_gc_balance):
                        best_tm_diff = tm_diff
                        best_gc_balance = gc_balance
                        best_primers = {
                            'mutation': f"{original_aa}{position}{desired_aa}",
                            'original_codon': original_codon,
                            'new_codon': desired_codon,
                            'forward_primer': forward_primer,
                            'reverse_primer': reverse_primer,
                            'forward_tm': round(fw_tm, 2),
                            'reverse_tm': round(rv_tm, 2),
                            'forward_gc': round(fw_gc, 2),
                            'reverse_gc': round(rv_gc, 2)
                        }

        if best_primers:
            results.append(best_primers)
        else:
            # If we couldn't find primers meeting all criteria, use the most conservative approach
            flank_length = 15  # Default flank length
            left_flank = dna_sequence[max(0, codon_pos - flank_length):codon_pos]
            right_flank = dna_sequence[codon_pos+3:min(len(dna_sequence), codon_pos + 3 + flank_length)]

            # Use the first codon option
            desired_codon = desired_codons[0] if desired_codons else 'NNN'

            # Create primers
            forward_primer = left_flank + desired_codon + right_flank
            reverse_primer = reverse_complement(forward_primer)

            fw_tm = mt.Tm_NN(Seq(forward_primer))
            rv_tm = mt.Tm_NN(Seq(reverse_primer))
            fw_gc = calculate_gc_content(forward_primer)
            rv_gc = calculate_gc_content(reverse_primer)
            
            results.append({
                'mutation': f"{original_aa}{position}{desired_aa}",
                'original_codon': original_codon,
                'new_codon': desired_codon,
                'forward_primer': forward_primer,
                'reverse_primer': reverse_primer,
                'forward_tm': round(fw_tm, 2),
                'reverse_tm': round(rv_tm, 2),
                'forward_gc': round(fw_gc, 2),
                'reverse_gc': round(rv_gc, 2),
                'warning': "Primer does not meet optimal criteria but should still work for mutagenesis."
            })

    return results

# Example usage
if __name__ == "__main__":
    # Sample GPCR DNA sequence (replace with your actual sequence)
    gpcr_sequence ="ATGAACGGGACCGCCAGCGTGGCGCTGTTCAACCTGGCCATTGCTGATCGCTACCTGGCCATCGTCCTCTCTGCCATCATGGGCAACATGCTGGTCATCAGCGCCATTGCCAACCCTATCATCTACTGCAAGGACCTTCACTACACTCTGATCACGCCCATTGTCATTGACGTGGTCACCTCAGTCGTCAACCCGCTCATCTACACTCTCTTCAGGACTTACGTGATCATGAGCATGGTGCCCTTTGGCCTGGTGGGTAATTCACTACTGGTCACTCCGTTCATCATGTGTCTTCAGATGAAACTGCCCGCCAAACGCCACCAAGGTATCAGCACAGAGACACGCCAGAACTTCCTCTCTTCATCATATCGCCTCACCTTGATGTTTAGGCTTAGGCAGTTCATCATTGCCTATGCCATTGTCGACTCCTTACTCGCCTTCTACCAGTCATTTGAAAATGTCATCTGCAAAGACCTGATTGTGTTCGTCTTTGGTCTCTGCATTGCCCTTTTCATCACGCCGCTCATCGTCATCGACAGGTACACATCCATTGCAAAAGCAGCACAATTTATCATCATCTGCTGGTTTACCATCCCATTCATCTACAGCCTCCGAAGCAGCTTCATATGCAGCTTAGCACTTGTGAATCCTTTCCACTTCCTCACTGGACAGTTGGCAGCCTGCAAGCAGATTTTCCACATCCTCAAACTGAAGTTACAAAGCAGTGATGCAAACGGTAAACAGCAGACCAAAGAGGCAAGCACCAGCCAGCAGGAAGATCCACAAGAACCTCTGGCTCCAGAACAGAGCGTCATCAAAGTATCCTTCGACTATTACTTTTTTCCCAAGAATACAGTTGAGTCAAAGTGTATCTTGTAG"

    # Define mutations as (original_aa, position, desired_aa)
    mutations = [
        ('F', 10, 'A'),  # F at position 10
        ('G', 86, 'D'),  # G at position 86
        ('Y', 138, 'S')  # Y at position 138
    ]

    # Design primers
    primer_results = design_mutagenesis_primers(gpcr_sequence, mutations)

    # Print results
    for result in primer_results:
        if 'error' in result:
            print(f"Mutation {result['mutation']}: {result['error']}")
        else:
            print(f"\nMutation: {result['mutation']}")
            if 'warning' in result:
                print(f"Warning: {result['warning']}")
            print(f"Original codon: {result['original_codon']} → New codon: {result['new_codon']}")
            print(f"Forward primer ({len(result['forward_primer'])} bp):")
            print(f"5'-{result['forward_primer']}-3'")
            print(f"Tm: {result['forward_tm']}°C, GC: {result['forward_gc']}%")
            print(f"Reverse primer ({len(result['reverse_primer'])} bp):")
            print(f"5'-{result['reverse_primer']}-3'")
            print(f"Tm: {result['reverse_tm']}°C, GC: {result['reverse_gc']}%")
