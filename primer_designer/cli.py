#!/usr/bin/env python
"""
Command-line interface for primer-designer
"""

import argparse
import json
import sys
from typing import List, Tuple

from primer_designer.core import design_primers
from primer_designer.sequence import validate_dna


def parse_mutations(mutation_strings: List[str]) -> List[Tuple[str, int, str]]:
    """Parse mutation strings into tuples"""
    mutations = []
    for mutation_str in mutation_strings:
        # Mutation format should be like "F10A" - original AA, position, new AA
        if len(mutation_str) < 3:
            print(f"Error: Invalid mutation format: {mutation_str}", file=sys.stderr)
            sys.exit(1)

        # Extract position (middle digits)
        pos_start = 1
        while pos_start < len(mutation_str) and not mutation_str[pos_start].isdigit():
            pos_start += 1

        pos_end = pos_start
        while pos_end < len(mutation_str) and mutation_str[pos_end].isdigit():
            pos_end += 1

        if pos_start == pos_end:
            print(f"Error: Could not parse position from mutation: {mutation_str}", file=sys.stderr)
            sys.exit(1)

        original_aa = mutation_str[:pos_start]
        position = int(mutation_str[pos_start:pos_end])
        desired_aa = mutation_str[pos_end:]

        if not original_aa or not desired_aa:
            print(f"Error: Invalid mutation format: {mutation_str}", file=sys.stderr)
            sys.exit(1)

        mutations.append((original_aa, position, desired_aa))

    return mutations


def main():
    """Main CLI entry point"""
    parser = argparse.ArgumentParser(description="Design primers for site-directed mutagenesis")
    parser.add_argument("sequence", help="DNA sequence to use for primer design")
    parser.add_argument("mutations", nargs="+", help="Mutations in format F10A (original AA, position, new AA)")
    parser.add_argument("--json", action="store_true", help="Output in JSON format")
    parser.add_argument("--min-flank-length", type=int, default=15, help="Minimum nucleotides flanking each side of mutation (default: 15)")
    parser.add_argument("--max-flank-length", type=int, default=18, help="Maximum nucleotides flanking each side of mutation (default: 18)")
    parser.add_argument("--min-tm", type=float, default=65.0, help="Minimum acceptable melting temperature in °C (default: 65.0)")
    parser.add_argument("--max-tm", type=float, default=75.0, help="Maximum acceptable melting temperature in °C (default: 75.0)")
    parser.add_argument("--optimal-tm", type=float, default=70.0, help="Target optimal melting temperature in °C (default: 70.0)")
    parser.add_argument("--max-tm-diff", type=float, default=5.0, help="Maximum allowed difference between forward/reverse Tm (default: 5.0)")
    parser.add_argument("--min-gc", type=float, default=35.0, help="Minimum acceptable GC content in percent (default: 35.0)")
    parser.add_argument("--max-gc", type=float, default=70.0, help="Maximum acceptable GC content in percent (default: 70.0)")
    parser.add_argument("--optimal-gc", type=float, default=50.0, help="Target optimal GC content in percent (default: 50.0)")

    args = parser.parse_args()

    # Parse mutations
    sequence = validate_dna(args.sequence)
    mutations = parse_mutations(args.mutations)

    # Design primers
    results = design_primers(
        sequence,
        mutations,
        min_flank_length=args.min_flank_length,
        max_flank_length=args.max_flank_length,
        min_tm=args.min_tm,
        max_tm=args.max_tm,
        optimal_tm=args.optimal_tm,
        max_tm_diff=args.max_tm_diff,
        min_gc=args.min_gc,
        max_gc=args.max_gc,
        optimal_gc=args.optimal_gc
    )

    # Output results
    if args.json:
        print(json.dumps(results, indent=2))
    else:
        for result in results:
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


if __name__ == "__main__":
    main()
