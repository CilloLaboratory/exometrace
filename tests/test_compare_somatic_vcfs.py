from pathlib import Path
import tempfile
import unittest

from scripts.compare_somatic_vcfs import load_variants


class SomaticComparisonTests(unittest.TestCase):
    def test_load_variants(self) -> None:
        variants = load_variants(Path("tests/test_data/vcf/deepsomatic_test.vcf.gz"), "deepsomatic")
        self.assertIn(("chr1", "5", "A", "T"), variants)
        self.assertEqual(variants[("chr1", "5", "A", "T")]["filter"], "PASS")


if __name__ == "__main__":
    unittest.main()
