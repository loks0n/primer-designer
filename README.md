# primer-designer

[![Tests](https://github.com/loks0n/primer-designer/actions/workflows/test.yml/badge.svg)](https://github.com/loks0n/primer-designer/actions/workflows/test.yml)
[![PyPI version](https://img.shields.io/pypi/v/primer-designer.svg)](https://pypi.org/project/primer-designer/)
[![Python Version](https://img.shields.io/pypi/pyversions/primer-designer.svg)](https://pypi.org/project/primer-designer/)
[![License](https://img.shields.io/github/license/loks0n/primer-designer.svg)](https://github.com/loks0n/primer-designer/blob/main/LICENSE)

A tool for designing site-directed mutagenesis primers for GPCR sequences or other protein-coding DNA sequences.

## Features

- Generates optimal mutagenesis primers based on standard molecular biology criteria
- Calculates melting temperatures using nearest-neighbor thermodynamic models
- Computes GC content for designed primers
- Validates mutation requests against the provided DNA sequence
- Selects the most suitable codon for the desired amino acid substitution
- Produces forward and reverse primers ready for ordering
- Optimizes parameters within configurable ranges (Tm, GC content, etc.)
- Provides fallback options for challenging sequences

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

- Python 3.13+
- FastAPI 0.115.12+
- uvicorn 0.34.0+
- typing-extensions 4.13.0+

## Usage

### As a Python Library

```python
from primer_designer.core import design_primers

# Your DNA sequence
dna_sequence = "ATGAACGGGACCGCCAGCGTGGCGCTGTTCAACCTGGCCATTGCTGATCGCTACCTGGCCATCGTCCTCTCTGCC..."

# Define mutations as (original_aa, position, desired_aa)
mutations = [
    ('F', 10, 'A'),  # Replace F at position 10 with A
    ('G', 86, 'D'),  # Replace G at position 86 with D
    ('Y', 138, 'S')  # Replace Y at position 138 with S
]

# Design primers with default parameters
primer_results = design_primers(dna_sequence, mutations)

# Design primers with custom parameters
primer_results = design_primers(
    dna_sequence, 
    mutations,
    min_flank_length=12,
    max_flank_length=20,
    min_tm=60.0,
    max_tm=80.0,
    optimal_tm=70.0,
    max_tm_diff=3.0,
    min_gc=40.0,
    max_gc=65.0,
    optimal_gc=55.0
)

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

# Custom parameters
primer-designer --min-tm 60.0 --max-gc 65.0 ATGAACGGGACC... F10A G86D
```

### REST API

The package includes a FastAPI server for providing primer design as a service:

```bash
# Start the API server
python -m primer_designer.api
```

API endpoint:
- POST `/primers` - Design primers based on the provided sequence and mutations

Example request:
```json
{
  "sequence": "ATGAACGGGACCGCCAGCGTGGCGCTGTTCAACCTGGCC...",
  "mutations": ["F10A", "G86D", "Y138S"],
  "min_flank_length": 15,
  "max_flank_length": 18,
  "min_tm": 65.0,
  "max_tm": 75.0,
  "optimal_tm": 70.0,
  "max_tm_diff": 5.0,
  "min_gc": 35.0,
  "max_gc": 70.0,
  "optimal_gc": 50.0
}
```

## How it works

Primer designer follows these steps to design optimal primers:

1. Translates the DNA sequence to verify the requested mutations
2. Identifies the codon position for each mutation
3. Selects the optimal codon for the desired amino acid
4. Designs primers with 15-18 nucleotides (configurable) flanking each side of the mutation
5. Calculates and optimizes for:
   - Melting temperature (Tm) ~70°C (configurable) with minimal difference between forward and reverse primers
   - GC content between 35-70% (configurable)
   - Overall primer length (typically 33-39 bp)
6. Falls back to relaxed criteria if optimal parameters can't be met

## Advanced Features

### Salt Correction for Melting Temperature

The melting temperature calculation includes sophisticated salt correction based on ion concentrations:
- Sodium (Na+), Potassium (K+), Tris buffer, and Magnesium (Mg2+) concentrations

### Nearest-Neighbor Thermodynamic Models

Primer melting temperatures are calculated using nearest-neighbor thermodynamic models (Allawi & SantaLucia, 1997), which provide more accurate Tm predictions than simpler %GC-based methods.

## Development

```bash
# Install development dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Run linting
ruff check .

# Run type checking
pyright
```

## License

[MIT License](LICENSE)