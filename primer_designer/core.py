from typing import List, Tuple

from primer_designer.sequence import (
    DNA,
    validate_dna,
    protein_to_codons,
    reverse_complement,
    gc_content,
    melting_temperature,
    translate
)

def design_primers(
    dna_sequence: DNA,
    mutations: List[Tuple[str, int, str]],
    min_flank_length: int = 15,
    max_flank_length: int = 18,
    min_tm: float = 65.0,
    max_tm: float = 75.0,
    optimal_tm: float = 70.0,
    max_tm_diff: float = 5.0,
    min_gc: float = 35.0,
    max_gc: float = 70.0,
    optimal_gc: float = 50.0
) -> List[dict]:
    """
    Design primers for site-directed mutagenesis with configurable criteria.

    Parameters:
    dna_sequence (str): The GPCR DNA sequence
    mutations (list): List of tuples (original_aa, position, desired_aa)
    min_flank_length (int): Minimum nucleotides flanking each side of mutation (default: 15)
    max_flank_length (int): Maximum nucleotides flanking each side of mutation (default: 18)
    min_tm (float): Minimum acceptable melting temperature in °C (default: 65.0)
    max_tm (float): Maximum acceptable melting temperature in °C (default: 75.0)
    optimal_tm (float): Target optimal melting temperature in °C (default: 70.0)
    max_tm_diff (float): Maximum allowed difference between forward/reverse Tm (default: 5.0)
    min_gc (float): Minimum acceptable GC content in percent (default: 35.0)
    max_gc (float): Maximum acceptable GC content in percent (default: 70.0)
    optimal_gc (float): Target optimal GC content in percent (default: 50.0)

    Returns:
    list: List of dictionaries containing mutation details and primers
    """
    protein_sequence = translate(dna_sequence)
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
        codon_pos = (position - 1) * 3
        original_codon = dna_sequence[codon_pos:codon_pos+3]

        # Generate all possible codons for the desired amino acid
        desired_codons = protein_to_codons(desired_aa)

        best_primers = None
        best_score = float('inf')

        for desired_codon in desired_codons:
            # Try different flank lengths from min to max
            for flank_length in range(min_flank_length, max_flank_length + 1):
                # Ensure we don't go out of bounds
                left_flank = dna_sequence[max(0, codon_pos - flank_length):codon_pos]
                right_flank = dna_sequence[codon_pos+3:min(len(dna_sequence), codon_pos + 3 + flank_length)]

                # Create primers with the mutation
                forward_primer = validate_dna(left_flank + desired_codon + right_flank)
                reverse_primer = reverse_complement(forward_primer)

                # Calculate properties
                fw_tm = melting_temperature(forward_primer)
                rv_tm = melting_temperature(reverse_primer)
                fw_gc = gc_content(forward_primer)
                rv_gc = gc_content(reverse_primer)

                # Check if primers meet requirements
                tm_requirement = min_tm <= fw_tm <= max_tm and min_tm <= rv_tm <= max_tm and abs(fw_tm - rv_tm) <= max_tm_diff
                gc_requirement = min_gc <= fw_gc <= max_gc and min_gc <= rv_gc <= max_gc

                if tm_requirement and gc_requirement:
                    # Calculate a combined score - lower is better
                    # Weight factors based on importance (can be adjusted)
                    tm_score = abs(fw_tm - optimal_tm) + abs(rv_tm - optimal_tm)
                    gc_score = abs(fw_gc - optimal_gc) + abs(rv_gc - optimal_gc)
                    total_score = tm_score * 1.5 + gc_score  # Tm is weighted more heavily

                    # Choose the best primers
                    if total_score < best_score:
                        best_score = total_score
                        best_primers = {
                            'mutation': f"{original_aa}{position}{desired_aa}",
                            'original_codon': original_codon,
                            'new_codon': desired_codon,
                            'forward_primer': forward_primer,
                            'reverse_primer': reverse_primer,
                            'forward_tm': round(fw_tm, 2),
                            'reverse_tm': round(rv_tm, 2),
                            'forward_gc': round(fw_gc, 3),
                            'reverse_gc': round(rv_gc, 3)
                        }

        if best_primers:
            results.append(best_primers)
        else:
            # If we couldn't find primers meeting all criteria, try with relaxed constraints
            relaxed_min_tm = min_tm - 5
            relaxed_max_tm = max_tm + 5
            relaxed_max_tm_diff = max_tm_diff + 3
            relaxed_min_gc = max(25, min_gc - 10)
            relaxed_max_gc = min(80, max_gc + 10)

            best_relaxed_primers = None
            best_relaxed_score = float('inf')

            for desired_codon in desired_codons:
                for flank_length in range(min_flank_length, max_flank_length + 1):
                    left_flank = dna_sequence[max(0, codon_pos - flank_length):codon_pos]
                    right_flank = dna_sequence[codon_pos+3:min(len(dna_sequence), codon_pos + 3 + flank_length)]

                    forward_primer = validate_dna(left_flank + desired_codon + right_flank)
                    reverse_primer = reverse_complement(forward_primer)

                    fw_tm = melting_temperature(forward_primer)
                    rv_tm = melting_temperature(reverse_primer)
                    fw_gc = gc_content(forward_primer)
                    rv_gc = gc_content(reverse_primer)

                    relaxed_tm_requirement = relaxed_min_tm <= fw_tm <= relaxed_max_tm and relaxed_min_tm <= rv_tm <= relaxed_max_tm and abs(fw_tm - rv_tm) <= relaxed_max_tm_diff
                    relaxed_gc_requirement = relaxed_min_gc <= fw_gc <= relaxed_max_gc and relaxed_min_gc <= rv_gc <= relaxed_max_gc

                    if relaxed_tm_requirement and relaxed_gc_requirement:
                        # Calculate a score for relaxed criteria
                        tm_score = abs(fw_tm - optimal_tm) + abs(rv_tm - optimal_tm)
                        gc_score = abs(fw_gc - optimal_gc) + abs(rv_gc - optimal_gc)
                        total_score = tm_score * 1.5 + gc_score

                        if total_score < best_relaxed_score:
                            best_relaxed_score = total_score
                            best_relaxed_primers = {
                                'mutation': f"{original_aa}{position}{desired_aa}",
                                'original_codon': original_codon,
                                'new_codon': desired_codon,
                                'forward_primer': forward_primer,
                                'reverse_primer': reverse_primer,
                                'forward_tm': round(fw_tm, 2),
                                'reverse_tm': round(rv_tm, 2),
                                'forward_gc': round(fw_gc, 3),
                                'reverse_gc': round(rv_gc, 3),
                                'warning': "Primer does not meet optimal criteria but may still work for mutagenesis."
                            }

            if best_relaxed_primers:
                results.append(best_relaxed_primers)
            else:
                # Last resort: create basic primers with default values
                flank_length = min_flank_length
                left_flank = dna_sequence[max(0, codon_pos - flank_length):codon_pos]
                right_flank = dna_sequence[codon_pos+3:min(len(dna_sequence), codon_pos + 3 + flank_length)]

                # Use the first codon option
                desired_codon = desired_codons[0] if desired_codons else 'NNN'

                # Create primers
                forward_primer = validate_dna(left_flank + desired_codon + right_flank)
                reverse_primer = reverse_complement(forward_primer)

                fw_tm = melting_temperature(forward_primer)
                rv_tm = melting_temperature(reverse_primer)
                fw_gc = gc_content(forward_primer)
                rv_gc = gc_content(reverse_primer)

                results.append({
                    'mutation': f"{original_aa}{position}{desired_aa}",
                    'original_codon': original_codon,
                    'new_codon': desired_codon,
                    'forward_primer': forward_primer,
                    'reverse_primer': reverse_primer,
                    'forward_tm': round(fw_tm, 2),
                    'reverse_tm': round(rv_tm, 2),
                    'forward_gc': round(fw_gc, 3),
                    'reverse_gc': round(rv_gc, 3),
                    'warning': "Primer design with fallback parameters - may not be optimal for mutagenesis."
                })

    return results

if __name__ == "__main__":
    # Sample GPCR DNA sequence (replace with your actual sequence)
    gpcr_sequence = validate_dna("ATGAACGGGACCGCCAGCGTGGCGCTGTTCAACCTGGCCATTGCTGATCGCTACCTGGCCATCGTCCTCTCTGCCATCATGGGCAACATGCTGGTCATCAGCGCCATTGCCAACCCTATCATCTACTGCAAGGACCTTCACTACACTCTGATCACGCCCATTGTCATTGACGTGGTCACCTCAGTCGTCAACCCGCTCATCTACACTCTCTTCAGGACTTACGTGATCATGAGCATGGTGCCCTTTGGCCTGGTGGGTAATTCACTACTGGTCACTCCGTTCATCATGTGTCTTCAGATGAAACTGCCCGCCAAACGCCACCAAGGTATCAGCACAGAGACACGCCAGAACTTCCTCTCTTCATCATATCGCCTCACCTTGATGTTTAGGCTTAGGCAGTTCATCATTGCCTATGCCATTGTCGACTCCTTACTCGCCTTCTACCAGTCATTTGAAAATGTCATCTGCAAAGACCTGATTGTGTTCGTCTTTGGTCTCTGCATTGCCCTTTTCATCACGCCGCTCATCGTCATCGACAGGTACACATCCATTGCAAAAGCAGCACAATTTATCATCATCTGCTGGTTTACCATCCCATTCATCTACAGCCTCCGAAGCAGCTTCATATGCAGCTTAGCACTTGTGAATCCTTTCCACTTCCTCACTGGACAGTTGGCAGCCTGCAAGCAGATTTTCCACATCCTCAAACTGAAGTTACAAAGCAGTGATGCAAACGGTAAACAGCAGACCAAAGAGGCAAGCACCAGCCAGCAGGAAGATCCACAAGAACCTCTGGCTCCAGAACAGAGCGTCATCAAAGTATCCTTCGACTATTACTTTTTTCCCAAGAATACAGTTGAGTCAAAGTGTATCTTGTAG")

    # Define mutations as (original_aa, position, desired_aa)
    mutations = [
        ('F', 10, 'A'),  # F at position 10
        ('G', 86, 'D'),  # G at position 86
        ('Y', 138, 'S')  # Y at position 138
    ]

    # Design primers
    primer_results = design_primers(gpcr_sequence, mutations)

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
            print(f"Tm: {result['forward_tm']}°C, GC: {result['forward_gc'] * 100}%")
            print(f"Reverse primer ({len(result['reverse_primer'])} bp):")
            print(f"5'-{result['reverse_primer']}-3'")
            print(f"Tm: {result['reverse_tm']}°C, GC: {result['reverse_gc'] * 100}%")
