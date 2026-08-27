from pathlib import Path
import tempfile
import unittest

from scripts.summarize_cnv import summarize, read_arms, read_segments


class CNVSummaryTests(unittest.TestCase):
    def test_arm_summary(self) -> None:
        rows = summarize(
            "P001",
            read_segments(Path("tests/test_data/cnv/example.call.cns")),
            read_arms(Path("tests/test_data/refs/GRCh38/annotations/chromosome_arms.bed")),
            {
                "deep_deletion_log2": -1.1,
                "loss_log2": -0.3,
                "gain_log2": 0.2,
                "amplification_log2": 0.7,
            },
        )
        self.assertEqual(rows[0]["status"], "loss")
        self.assertEqual(rows[1]["status"], "gain")

    def test_read_segments_allows_large_fields(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "large.call.cns"
            large_gene_field = "GENE" * 40000
            path.write_text(
                "chromosome\tstart\tend\tgene\tlog2\n"
                f"chr1\t0\t100\t{large_gene_field}\t0.25\n",
                encoding="utf-8",
            )

            rows = read_segments(path)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["chromosome"], "chr1")
        self.assertEqual(rows[0]["log2"], "0.25")


if __name__ == "__main__":
    unittest.main()
