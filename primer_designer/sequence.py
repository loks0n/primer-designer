from typing import Dict, Tuple
from typing import NewType
import math
import re

DNA = NewType("DNA", str)
Proteins = NewType("Proteins", str)

# Thermodynamic lookup tables for nearest neighbors
# Values are (enthalpy, entropy) in (kcal/mol, cal/mol K)
# Default table: Allawi & SantaLucia (1997)
DNA_NN3 = {
    "init": (0, 0),
    "init_A/T": (2.3, 4.1),
    "init_G/C": (0.1, -2.8),
    "init_oneG/C": (0, 0),
    "init_allA/T": (0, 0),
    "init_5T/A": (0, 0),
    "sym": (0, -1.4),
    "AA/TT": (-7.9, -22.2),
    "AT/TA": (-7.2, -20.4),
    "TA/AT": (-7.2, -21.3),
    "CA/GT": (-8.5, -22.7),
    "GT/CA": (-8.4, -22.4),
    "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2),
    "CG/GC": (-10.6, -27.2),
    "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.9),
}


def protein_to_codons(amino_acid):
    """Get all possible codons for a given amino acid."""
    amino_acid = amino_acid.upper()
    codon_map = {
        "A": ["GCT", "GCC", "GCA", "GCG"],
        "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "N": ["AAT", "AAC"],
        "D": ["GAT", "GAC"],
        "C": ["TGT", "TGC"],
        "Q": ["CAA", "CAG"],
        "E": ["GAA", "GAG"],
        "G": ["GGT", "GGC", "GGA", "GGG"],
        "H": ["CAT", "CAC"],
        "I": ["ATT", "ATC", "ATA"],
        "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
        "K": ["AAA", "AAG"],
        "M": ["ATG"],
        "F": ["TTT", "TTC"],
        "P": ["CCT", "CCC", "CCA", "CCG"],
        "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
        "T": ["ACT", "ACC", "ACA", "ACG"],
        "W": ["TGG"],
        "Y": ["TAT", "TAC"],
        "V": ["GTT", "GTC", "GTA", "GTG"],
        "*": ["TAA", "TAG", "TGA"],
    }
    if amino_acid not in codon_map:
        raise ValueError(f"Invalid amino acid: {amino_acid}")
    return codon_map[amino_acid]


def validate_dna(sequence: str) -> DNA:
    """
    Validate that the input string is a valid DNA sequence.

    Args:
        sequence (str): Input sequence to validate

    Returns:
        DNA: Validated DNA sequence (uppercase)

    Raises:
        ValueError: If sequence contains invalid characters
    """
    # Convert to uppercase and remove whitespace
    cleaned_seq = sequence.upper().replace(" ", "")

    # Check if sequence contains only valid DNA bases
    if not re.match(r"^[ACGT]+$", cleaned_seq):
        raise ValueError("Invalid DNA sequence. Only A, C, G, T are allowed.")

    return DNA(cleaned_seq)


def gc_content(dna_seq: DNA) -> float:
    """
    Calculate GC content of a DNA sequence.

    Args:
        dna_seq (DNA): Input DNA sequence

    Returns:
        float: Percentage of G and C bases (0.0 to 1.0)
    """
    if not dna_seq:
        return 0.0

    gc_count = dna_seq.count("G") + dna_seq.count("C")
    return gc_count / len(dna_seq)


def reverse_complement(dna_seq: DNA) -> DNA:
    """
    Reverse complement a DNA sequence.

    Args:
        dna_seq (DNA): Input DNA sequence

    Returns:
        DNA: Reverse complement of the input sequence
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return DNA("".join(complement[base] for base in reversed(dna_seq)))


def translate(dna_seq: DNA) -> Proteins:
    codon_to_protein_dict = {
        "ATA": "I",
        "ATC": "I",
        "ATT": "I",
        "ATG": "M",
        "ACA": "T",
        "ACC": "T",
        "ACG": "T",
        "ACT": "T",
        "AAC": "N",
        "AAT": "N",
        "AAA": "K",
        "AAG": "K",
        "AGC": "S",
        "AGT": "S",
        "AGA": "R",
        "AGG": "R",
        "CTA": "L",
        "CTC": "L",
        "CTG": "L",
        "CTT": "L",
        "CCA": "P",
        "CCC": "P",
        "CCG": "P",
        "CCT": "P",
        "CAC": "H",
        "CAT": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGA": "R",
        "CGC": "R",
        "CGG": "R",
        "CGT": "R",
        "GTA": "V",
        "GTC": "V",
        "GTG": "V",
        "GTT": "V",
        "GCA": "A",
        "GCC": "A",
        "GCG": "A",
        "GCT": "A",
        "GAC": "D",
        "GAT": "D",
        "GAA": "E",
        "GAG": "E",
        "GGA": "G",
        "GGC": "G",
        "GGG": "G",
        "GGT": "G",
        "TCA": "S",
        "TCC": "S",
        "TCG": "S",
        "TCT": "S",
        "TTC": "F",
        "TTT": "F",
        "TTA": "L",
        "TTG": "L",
        "TAC": "Y",
        "TAT": "Y",
        "TAA": "*",
        "TAG": "*",
        "TGC": "C",
        "TGT": "C",
        "TGA": "*",
        "TGG": "W",
    }
    proteins = ""
    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i : i + 3]
        if len(codon) != 3 or codon not in codon_to_protein_dict:
            raise ValueError(f"Invalid codon: {codon}")
        proteins += codon_to_protein_dict[codon]
    return Proteins(proteins)


def salt_correction(
    seq: DNA, Na: float = 50, K: float = 0, Tris: float = 0, Mg: float = 0
) -> float:
    """
    Calculate salt correction term for melting temperature using method 5.

    Args:
        Na: Sodium concentration [mM]
        K: Potassium concentration [mM]
        Tris: Tris buffer concentration [mM]
        Mg: Magnesium concentration [mM]
        seq: DNA sequence

    Returns:
        float: Salt correction term
    """
    if not seq:
        raise ValueError("Sequence is required for salt correction method 5")

    # Validate inputs
    if Na < 0 or K < 0 or Tris < 0 or Mg < 0:
        raise ValueError("Ion concentrations cannot be negative")

    # Calculate monovalent ion concentration
    Mon = Na + K + Tris / 2.0  # millimolar concentration
    mon = Mon * 1e-3  # Convert to molar

    # Na equivalent according to von Ahsen et al. (2001)
    if Mg > 0:
        # Apply the magnesium correction to mon directly
        # to ensure it affects the final correction value
        mon += 120 * math.sqrt(Mg * 1e-3) * 1e-3

    if not mon:
        raise ValueError("Total ion concentration of zero is not allowed")

    # Correction
    return 0.368 * (len(seq) - 1) * math.log(mon)


def melting_temperature(
    dna_seq: DNA,
    nn_table: Dict[str, Tuple[float, float]] = DNA_NN3,
    Na: float = 50,
    K: float = 0,
    Tris: float = 0,
    Mg: float = 0,
    dnac1: float = 25,
    dnac2: float = 25,
    selfcomp: bool = False,
) -> float:
    """
    Calculate the melting temperature using the nearest neighbor method with salt correction.

    Args:
        dna_seq (DNA): Input DNA sequence
        nn_table: Thermodynamic nearest neighbor table
        Na: Sodium concentration [mM]
        K: Potassium concentration [mM]
        Tris: Tris buffer concentration [mM]
        Mg: Magnesium concentration [mM]
        dnac1: Concentration of the higher concentrated strand [nM]
        dnac2: Concentration of the lower concentrated strand [nM]
        selfcomp: Is the sequence self-complementary?

    Returns:
        float: Melting temperature in Celsius
    """

    # Create complementary sequence
    comp_map = {"A": "T", "T": "A", "G": "C", "C": "G"}
    c_seq = "".join(comp_map[base] for base in dna_seq)

    delta_h = 0  # Enthalpy
    delta_s = 0  # Entropy

    # Indexes for enthalpy and entropy in the thermodynamic tables
    d_h = 0
    d_s = 1

    # General initiation value
    delta_h += nn_table["init"][d_h]
    delta_s += nn_table["init"][d_s]

    # Check if sequence has all A/T or at least one G/C
    if gc_content(dna_seq) == 0:
        delta_h += nn_table["init_allA/T"][d_h]
        delta_s += nn_table["init_allA/T"][d_s]
    else:
        delta_h += nn_table["init_oneG/C"][d_h]
        delta_s += nn_table["init_oneG/C"][d_s]

    # Penalty if 5' end is T
    if dna_seq.startswith("T"):
        delta_h += nn_table["init_5T/A"][d_h]
        delta_s += nn_table["init_5T/A"][d_s]
    if dna_seq.endswith("A"):
        delta_h += nn_table["init_5T/A"][d_h]
        delta_s += nn_table["init_5T/A"][d_s]

    # Terminal basepairs
    ends = dna_seq[0] + dna_seq[-1]
    AT = ends.count("A") + ends.count("T")
    GC = ends.count("G") + ends.count("C")
    delta_h += nn_table["init_A/T"][d_h] * AT
    delta_s += nn_table["init_A/T"][d_s] * AT
    delta_h += nn_table["init_G/C"][d_h] * GC
    delta_s += nn_table["init_G/C"][d_s] * GC

    # Calculate nearest neighbor thermodynamics
    for i in range(len(dna_seq) - 1):
        nn_pair = dna_seq[i : i + 2] + "/" + c_seq[i : i + 2]

        # Check if the pair exists in the table
        if nn_pair in nn_table:
            delta_h += nn_table[nn_pair][d_h]
            delta_s += nn_table[nn_pair][d_s]
        elif nn_pair[::-1] in nn_table:  # Try reverse
            delta_h += nn_table[nn_pair[::-1]][d_h]
            delta_s += nn_table[nn_pair[::-1]][d_s]
        else:
            # If pair not found, try using symmetry rule
            reversed_pair = nn_pair[3:] + nn_pair[2] + "/" + nn_pair[0:2]
            if reversed_pair in nn_table:
                delta_h += nn_table[reversed_pair][d_h]
                delta_s += nn_table[reversed_pair][d_s]
            else:
                # If still not found, use approximation
                print(
                    f"Warning: thermodynamic data for {nn_pair} not found. Using approximation."
                )
                # Use average values as approximation
                delta_h += -8.0
                delta_s += -22.0

    # Apply salt correction
    salt_corr = salt_correction(dna_seq, Na=Na, K=K, Tris=Tris, Mg=Mg)
    delta_s += salt_corr

    # Calculate the concentration term
    R = 1.987  # Universal gas constant in Cal/K*mol
    k = (dnac1 - (dnac2 / 2.0)) * 1e-9
    if selfcomp:
        k = dnac1 * 1e-9
        delta_h += nn_table["sym"][d_h]
        delta_s += nn_table["sym"][d_s]

    # Calculate melting temperature
    melting_temp = (1000 * delta_h) / (delta_s + (R * math.log(k))) - 273.15

    return melting_temp
