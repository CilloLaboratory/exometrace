from pathlib import Path
import tempfile
import unittest

from scripts.merge_sample_qc import main as merge_main
from scripts.summarize_alignment_qc import parse_duplication_metrics, parse_flagstat, parse_hs_metrics, parse_insert_metrics
from scripts.validate_fastq_pairs import validate_pair


class FastqValidationTests(unittest.TestCase):
    def test_validate_pair(self) -> None:
        r1 = Path("tests/test_data/fastq/P001_T_R1.fastq.gz")
        r2 = Path("tests/test_data/fastq/P001_T_R2.fastq.gz")
        reads_r1, reads_r2 = validate_pair(r1, r2)
        self.assertEqual((reads_r1, reads_r2), (2, 2))


class QCSummaryParsersTests(unittest.TestCase):
    def test_flagstat_parser(self) -> None:
        path = Path("tests/test_data/qc/example.flagstat.txt")
        total, mapped, mapping_rate, proper_pair_rate = parse_flagstat(path)
        self.assertEqual(total, 100)
        self.assertEqual(mapped, 98)
        self.assertAlmostEqual(mapping_rate, 0.98)
        self.assertAlmostEqual(proper_pair_rate, 0.95)

    def test_picard_and_hs_parsers(self) -> None:
        dup_path = Path("tests/test_data/qc/example.duplication_metrics.txt")
        ins_path = Path("tests/test_data/qc/example.insert_size_metrics.txt")
        hs_path = Path("tests/test_data/qc/example.hs_metrics.txt")
        self.assertAlmostEqual(parse_duplication_metrics(dup_path), 0.1)
        self.assertAlmostEqual(parse_insert_metrics(ins_path), 250.0)
        hs = parse_hs_metrics(hs_path)
        self.assertAlmostEqual(hs["mean_target_coverage"], 120.0)
        self.assertAlmostEqual(hs["pct_target_20x"], 0.97)


if __name__ == "__main__":
    unittest.main()
