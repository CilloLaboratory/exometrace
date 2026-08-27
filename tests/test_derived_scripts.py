import unittest

from pathlib import Path
import tempfile

from scripts.build_callable_bed import depth_to_intervals, total_bases
from scripts.calculate_ccf import classify_clonality, expected_mutant_copies, find_segment, infer_ccf
from scripts.calculate_signature_exposure import (
    parse_activity_table,
    resolve_sigprofiler_volume,
    summarize_exposures,
    write_vcf_from_maf,
)
from scripts.calculate_tmb import calculate_tmb


class DerivedMetricsTests(unittest.TestCase):
    def test_callable_bases(self) -> None:
        bed = Path("tests/test_data/refs/GRCh38/intervals/exome_targets.bed")
        self.assertEqual(total_bases(bed), 12)

    def test_depth_to_intervals(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            depth_path = Path(tmpdir) / "callable.depth.tsv"
            depth_path.write_text(
                "\n".join(
                    [
                        "chr1\t1\t25\t15",
                        "chr1\t2\t25\t15",
                        "chr1\t3\t10\t15",
                        "chr1\t4\t25\t15",
                        "chr1\t5\t25\t5",
                        "chr1\t6\t25\t15",
                        "chr2\t1\t30\t12",
                        "chr2\t2\t30\t12",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            self.assertEqual(
                depth_to_intervals(depth_path, tumor_min_depth=20, normal_min_depth=10),
                [("chr1", 0, 2), ("chr1", 3, 4), ("chr1", 5, 6), ("chr2", 0, 2)],
            )

    def test_tmb_calculation(self) -> None:
        self.assertAlmostEqual(calculate_tmb(2, 0.000012), 166666.66666666666)

    def test_sigprofiler_activity_parser_row_oriented(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            table = Path(tmpdir) / "activities.tsv"
            table.write_text(
                "Samples\tSBS1\tSBS2\tSBS7a\tSBS13\n"
                "P001\t12\t3\t4\t2\n",
                encoding="utf-8",
            )
            self.assertEqual(
                parse_activity_table(table, "P001"),
                {"SBS1": 12.0, "SBS2": 3.0, "SBS7a": 4.0, "SBS13": 2.0},
            )

    def test_sigprofiler_activity_parser_column_oriented(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            table = Path(tmpdir) / "activities.tsv"
            table.write_text(
                "Signature\tP001\n"
                "SBS1\t12\n"
                "SBS7a\t4\n",
                encoding="utf-8",
            )
            self.assertEqual(
                parse_activity_table(table, "P001"),
                {"SBS1": 12.0, "SBS7a": 4.0},
            )

    def test_signature_summary_collapses_uv_and_apobec(self) -> None:
        exposure_row, qc_row = summarize_exposures(
            "P001",
            75,
            {"SBS1": 10.0, "SBS2": 3.0, "SBS7a": 4.0, "SBS13": 2.0},
            min_mutations_for_qc_pass=50,
        )
        self.assertEqual(exposure_row["SBS7_UV"], "4.000000")
        self.assertEqual(exposure_row["SBS2_APOBEC"], "5.000000")
        self.assertEqual(qc_row["status"], "PASS")

    def test_write_vcf_from_maf_uses_requested_contig_style(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_vcf = Path(tmpdir) / "P001.vcf"
            write_vcf_from_maf(
                [
                    {
                        "Chromosome": "1",
                        "Start_Position": "5",
                        "Reference_Allele": "A",
                        "Tumor_Seq_Allele2": "T",
                    }
                ],
                output_vcf,
                "P001",
                "chr",
            )
            self.assertIn("chr1\t5\t.\tA\tT", output_vcf.read_text(encoding="utf-8"))

    def test_resolve_sigprofiler_volume_uses_reference_config_base_dir(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            config_dir = root / "config"
            config_dir.mkdir(parents=True, exist_ok=True)
            reference_config = {"sigprofiler": {"volume_dir": "references/GRCh38/sigprofiler"}}
            pipeline_config = {"signatures": {"volume": None}}
            sig_dir = root / "references" / "GRCh38" / "sigprofiler" / "tsb" / "GRCh38"
            sig_dir.mkdir(parents=True, exist_ok=True)
            reference_config_path = config_dir / "references.yaml"
            pipeline_config_path = config_dir / "default.yaml"
            self.assertEqual(
                resolve_sigprofiler_volume(
                    pipeline_config,
                    pipeline_config_path,
                    reference_config,
                    reference_config_path,
                ),
                str((root / "references" / "GRCh38" / "sigprofiler").resolve()),
            )

    def test_resolve_sigprofiler_volume_rejects_missing_assets(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            config_dir = root / "config"
            config_dir.mkdir(parents=True, exist_ok=True)
            reference_config = {"sigprofiler": {"volume_dir": "references/GRCh38/sigprofiler"}}
            pipeline_config = {"signatures": {"volume": None}}
            sig_dir = root / "references" / "GRCh38" / "sigprofiler"
            sig_dir.mkdir(parents=True, exist_ok=True)
            with self.assertRaises(FileNotFoundError):
                resolve_sigprofiler_volume(
                    pipeline_config,
                    config_dir / "default.yaml",
                    reference_config,
                    config_dir / "references.yaml",
                )

    def test_find_segment_matches_chr_prefix_independently(self) -> None:
        segment = find_segment(
            [
                {"chromosome": "chr1", "start": "0", "end": "10", "total_cn": "2", "major_cn": "1", "minor_cn": "1"},
            ],
            "1",
            5,
        )
        self.assertIsNotNone(segment)
        self.assertEqual(segment["total_cn"], "2")

    def test_expected_mutant_copies_adjusts_for_purity_and_copy_number(self) -> None:
        self.assertEqual(expected_mutant_copies(vaf=0.25, total_cn=4.0, purity=0.75), 1)

    def test_infer_ccf_supports_clonal_interpretation(self) -> None:
        ccf, ccf_lower, ccf_upper, prob_clonal = infer_ccf(
            t_alt_count=30,
            t_depth=40,
            purity=0.75,
            total_cn=2.0,
            mutant_copies=1,
        )
        self.assertGreaterEqual(ccf, 0.9)
        self.assertGreater(prob_clonal, 0.5)
        self.assertEqual(classify_clonality(ccf, ccf_lower, ccf_upper, 0.9), "clonal")

    def test_infer_ccf_supports_subclonal_interpretation(self) -> None:
        ccf, ccf_lower, ccf_upper, _prob_clonal = infer_ccf(
            t_alt_count=8,
            t_depth=60,
            purity=0.75,
            total_cn=2.0,
            mutant_copies=1,
        )
        self.assertLess(ccf_upper, 0.9)
        self.assertEqual(classify_clonality(ccf, ccf_lower, ccf_upper, 0.9), "subclonal")


if __name__ == "__main__":
    unittest.main()
