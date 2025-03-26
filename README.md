# primer-designer

[![Tests](https://github.com/loks0n/primer-designer/actions/workflows/test.yml/badge.svg)](https://github.com/loks0n/primer-designer/actions/workflows/test.yml)
[![PyPI version](https://img.shields.io/pypi/v/primer-designer.svg)](https://pypi.org/project/primer-designer/)
[![Python Version](https://img.shields.io/pypi/pyversions/primer-designer.svg)](https://pypi.org/project/primer-designer/)
[![License](https://img.shields.io/github/license/loks0n/primer-designer.svg)](https://github.com/loks0n/primer-designer/blob/main/LICENSE)

A tool for designing site-directed mutagenesis primers for GPCR sequences or other protein-coding DNA sequences.

## Features

- Generates optimal mutagenesis primers based on standard molecular biology criteria
- Calculates melting temperatures and GC content for designed primers
- Validates mutation requests against the provided DNA sequence
- Selects the most suitable codon for the desired amino acid substitution
- Produces forward and reverse primers ready for ordering

## Installation

### From PyPI

```bash
pip install primer-designer
```

### From Source

```bash
# Clone the repository
git clone https://github.com/loks0n/primer-designer.git
cd primer-designer

# Install in development mode
pip install -e .
```

## Requirements

- Python 3.13
- Biopython 1.85

## Usage

### As a Python Library

```python
from primer_designer.core import design_mutagenesis_primers

# Your DNA sequence
dna_sequence = "ATGAACGGGACCGCCAGCGTGGCGCTGTTCAACCTGGCCATTGCTGATCGCTACCTGGCCATCGTCCTCTCTGCC..."

# Define mutations as (original_aa, position, desired_aa)
mutations = [
    ('F', 10, 'A'),  # Replace F at position 10 with A
    ('G', 86, 'D'),  # Replace G at position 86 with D
    ('Y', 138, 'S')  # Replace Y at position 138 with S
]

# Design primers
primer_results = design_mutagenesis_primers(dna_sequence, mutations)

# Print results
for result in primer_results:
    print(f"\nMutation: {result['mutation']}")
    print(f"Forward primer: 5'-{result['forward_primer']}-3'")
    print(f"Reverse primer: 5'-{result['reverse_primer']}-3'")
```

### Command Line Interface

```bash
# Basic usage
primer-designer ATGAACGGGACC... F10A G86D Y138S

# Output as JSON
primer-designer --json ATGAACGGGACC... F10A G86D Y138S
```

## How it works

Primer designer follows these steps to design optimal primers:

1. Translates the DNA sequence to verify the requested mutations
2. Identifies the codon position for each mutation
3. Selects the optimal codon for the desired amino acid
4. Designs primers with 15-18 nucleotides flanking each side of the mutation
5. Calculates and optimizes for:
   - Melting temperature (Tm) ~65°C with minimal difference between forward and reverse primers
   - GC content between 35-70%
   - Overall primer length (typically 33-39 bp)

## Example Output

```
Mutation: F10A
Original codon: TTC → New codon: GCC
Forward primer (33 bp):
5'-ATGAACGGGACCGCCAGCGTGGCGCTGTTCGCCCTGGC-3'
Tm: 72.11°C, GC: 68.97%
Reverse primer (33 bp):
5'-GCCAGGGCGAACAGCGCCACGCTGGCGGTCCCGTTCAT-3'
Tm: 72.11°C, GC: 68.97%
```

## License

[MIT License](LICENSE)