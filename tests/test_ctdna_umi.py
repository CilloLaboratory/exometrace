import csv
import gzip
import tempfile
import unittest

from pathlib import Path

from scripts.trim_ctdna_template import main as trim_ctdna_template_main
from scripts.trim_ctdna_template import parse_read_structure


REPO_ROOT = Path(__file__).resolve().parents[1]


class CtDnaTemplateTrimTests(unittest.TestCase):
    def test_parse_read_structure_returns_roche_trim_length(self) -> None:
        parsed = parse_read_structure("3M3S+T")
        self.assertEqual(parsed.template_trim_bases, 6)

    def test_parse_read_structure_rejects_unsupported_shapes(self) -> None:
        for value in ("", "8B+T", "3M+T2S", "T", "3M3S"):
            with self.assertRaises(ValueError):
                parse_read_structure(value)

    def test_template_trim_script_removes_non_template_cycles(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            temp_root = Path(tempdir)
            r1 = temp_root / "input_R1.fastq.gz"
            r2 = temp_root / "input_R2.fastq.gz"
            output_r1 = temp_root / "trimmed_R1.fastq.gz"
            output_r2 = temp_root / "trimmed_R2.fastq.gz"
            qc_output = temp_root / "trim.tsv"

            with gzip.open(r1, "wt", encoding="utf-8") as handle:
                handle.write("@read1\nAAACCCGGTT\n+\nJJJJJJJJJJ\n")
            with gzip.open(r2, "wt", encoding="utf-8") as handle:
                handle.write("@read1\nTTTGGGCCAA\n+\nHHHHHHHHHH\n")

            import sys

            argv = sys.argv
            sys.argv = [
                "trim_ctdna_template.py",
                "--r1", str(r1),
                "--r2", str(r2),
                "--read-structure-r1", "3M3S+T",
                "--read-structure-r2", "3M3S+T",
                "--output-r1", str(output_r1),
                "--output-r2", str(output_r2),
                "--qc-output", str(qc_output),
                "--sample-id", "PLASMA",
            ]
            try:
                self.assertEqual(trim_ctdna_template_main(), 0)
            finally:
                sys.argv = argv

            with gzip.open(output_r1, "rt", encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "@read1\nGGTT\n+\nJJJJ\n")
            with gzip.open(output_r2, "rt", encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "@read1\nCCAA\n+\nHHHH\n")

            with qc_output.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["read_structure_r1"], "3M3S+T")
            self.assertEqual(rows[0]["template_trim_bases_r1"], "6")
            self.assertEqual(rows[0]["template_trim_bases_r2"], "6")


class CtDnaWorkflowShapeTests(unittest.TestCase):
    def test_defaults_use_roche_read_structures(self) -> None:
        config_text = (REPO_ROOT / "config" / "default.yaml").read_text(encoding="utf-8")
        self.assertIn("read_structure_r1: 3M3S+T", config_text)
        self.assertIn("read_structure_r2: 3M3S+T", config_text)
        self.assertNotIn("inline_r1_bases", config_text)
        self.assertNotIn("inline_r2_bases", config_text)

    def test_consensus_module_extracts_umis_with_both_read_structures(self) -> None:
        module_text = (REPO_ROOT / "modules" / "umi_consensus" / "main.nf").read_text(encoding="utf-8")
        self.assertIn("fgbio ExtractUmisFromBam", module_text)
        self.assertIn("-r ${read_structure_r1} ${read_structure_r2}", module_text)
        self.assertIn("read_structure_r1", module_text)
        self.assertIn("read_structure_r2", module_text)

    def test_main_routes_trimmed_fastqs_only_to_cfsnv_path(self) -> None:
        workflow_text = (REPO_ROOT / "main.nf").read_text(encoding="utf-8")
        self.assertIn("template_trimmed = UMI_TEMPLATE_TRIM", workflow_text)
        self.assertIn("ctdna_align_inputs = fastqc_results.map", workflow_text)
        self.assertIn("consensus_inputs = ctdna_aligned.map", workflow_text)
        self.assertIn("ctdna_read_structure_r1", workflow_text)
        self.assertIn("ctdna_read_structure_r2", workflow_text)


if __name__ == "__main__":
    unittest.main()
