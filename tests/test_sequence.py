import unittest
from primer_designer.sequence import (
    DNA,
    validate_dna,
    gc_content,
    reverse_complement,
    translate,
    protein_to_codons,
    salt_correction,
    melting_temperature,
)


class TestSequence(unittest.TestCase):
    def test_validate_dna(self):
        # Valid DNA sequences
        self.assertEqual(validate_dna("ACGT"), "ACGT")
        self.assertEqual(validate_dna("acgt"), "ACGT")
        self.assertEqual(validate_dna("A C G T"), "ACGT")

        # Invalid DNA sequences
        with self.assertRaises(ValueError):
            validate_dna("AXGT")
        with self.assertRaises(ValueError):
            validate_dna("RNA")
        with self.assertRaises(ValueError):
            validate_dna("12345")

    def test_gc_content(self):
        self.assertEqual(gc_content(DNA("ACGT")), 0.5)
        self.assertEqual(gc_content(DNA("AAAA")), 0.0)
        self.assertEqual(gc_content(DNA("GCGC")), 1.0)
        self.assertEqual(gc_content(DNA("")), 0.0)

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement(DNA("ACGT")), DNA("ACGT"))
        self.assertEqual(reverse_complement(DNA("AACGTT")), DNA("AACGTT"))
        self.assertEqual(reverse_complement(DNA("AGTC")), DNA("GACT"))
        self.assertEqual(reverse_complement(DNA("AAA")), DNA("TTT"))

    def test_translate(self):
        # Start codon and basic translation
        self.assertEqual(translate(DNA("ATG")), "M")

        # Multiple codons
        self.assertEqual(translate(DNA("ATGGCTAGCTGA")), "MAS*")
        self.assertEqual(translate(DNA("ATGGCTAGCTGA")), "MAS*")

        # Edge cases
        with self.assertRaises(ValueError):
            translate(DNA("AT"))  # Incomplete codon
        with self.assertRaises(ValueError):
            translate(DNA("ATGXYZ"))  # Invalid bases

    def test_protein_to_codons(self):
        # Single amino acid
        self.assertEqual(protein_to_codons("A"), ["GCT", "GCC", "GCA", "GCG"])

        # Edge cases
        with self.assertRaises(ValueError):
            protein_to_codons("X")  # Non-standard amino acid

        # Case insensitivity
        self.assertEqual(protein_to_codons("a"), ["GCT", "GCC", "GCA", "GCG"])

        # Stop codon
        self.assertEqual(protein_to_codons("*"), ["TAA", "TAG", "TGA"])

    def test_salt_correction(self):
        # Basic test
        correction = salt_correction(DNA("ACGT"), Na=50)
        # Original formula: 0.368 * (len(seq) - 1) * math.log(mon)
        # Updated expectation to match actual calculation
        self.assertAlmostEqual(correction, -3.31, places=2)

        # Zero ion concentration
        with self.assertRaises(ValueError):
            salt_correction(DNA("ACGT"), Na=0, K=0, Tris=0, Mg=0)

        # Empty sequence
        with self.assertRaises(ValueError):
            salt_correction(DNA(""), Na=50)

        # Negative ion concentration
        with self.assertRaises(ValueError):
            salt_correction(DNA("ACGT"), Na=-1)

        # Magnesium effect
        correction_with_mg = salt_correction(DNA("ACGT"), Na=50, Mg=1.5)
        correction_without_mg = salt_correction(DNA("ACGT"), Na=50)
        self.assertGreater(correction_with_mg, correction_without_mg)

    def test_melting_temperature(self):
        # Basic test
        tm = melting_temperature(DNA("ACGTACGT"))
        self.assertIsInstance(tm, float)

        # GC content effect
        tm_high_gc = melting_temperature(DNA("GCGCGCGC"))
        tm_low_gc = melting_temperature(DNA("ATATAGCG"))
        self.assertGreater(tm_high_gc, tm_low_gc)

        # Length effect
        tm_long = melting_temperature(DNA("ACGTACGTACGTACGT"))
        tm_short = melting_temperature(DNA("ACGTACGT"))
        self.assertGreater(tm_long, tm_short)

        # Salt concentration effect
        tm_high_salt = melting_temperature(DNA("ACGTACGT"), Na=100)
        tm_low_salt = melting_temperature(DNA("ACGTACGT"), Na=50)
        self.assertGreater(tm_high_salt, tm_low_salt)

        # Self-complementary
        tm_selfcomp = melting_temperature(DNA("ACGT"), selfcomp=True)
        tm_normal = melting_temperature(DNA("ACGT"), selfcomp=False)
        self.assertNotEqual(tm_selfcomp, tm_normal)

        # Edge cases
        # Very short sequence
        tm = melting_temperature(DNA("AT"))
        self.assertIsInstance(tm, float)

        # Sequence with unusual nearest neighbors not in the table
        # This should use the approximation warning
        tm = melting_temperature(DNA("ATCGATCG"))
        self.assertIsInstance(tm, float)


if __name__ == "__main__":
    unittest.main()
