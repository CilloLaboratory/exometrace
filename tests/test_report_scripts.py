from pathlib import Path
import csv
import tempfile
import unittest

from scripts.generate_containers_manifest import main as containers_main
from scripts.generate_execution_report import main as execution_report_main
from scripts.generate_multiqc_stub import main as multiqc_main


class ReportScriptTests(unittest.TestCase):
    def test_container_manifest_generation(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            output = Path(tempdir) / "containers.tsv"
            import sys
            argv = sys.argv
            sys.argv = [
                "generate_containers_manifest.py",
                "--config", "config/containers.yaml",
                "--output", str(output),
            ]
            try:
                containers_main()
            finally:
                sys.argv = argv
            with output.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["tool"] == "deepsomatic_gpu" for row in rows))

    def test_multiqc_dashboard_generation(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            sample_qc = Path(tempdir) / "sample_qc.tsv"
            output = Path(tempdir) / "multiqc_report.html"
            sample_qc.write_text(
                "patient_id\tsample\tsample_type\ttotal_reads\tmapped_reads\tmapping_rate\tproper_pair_rate\tduplicate_rate\tmean_target_coverage\tmedian_target_coverage\tpct_target_10x\tpct_target_20x\tpct_target_30x\tpct_target_50x\tpct_target_100x\tmean_insert_size\n"
                "P001\tP001_T\ttumor\t100\t98\t0.98\t0.95\t0.10\t120\t110\t0.99\t0.97\t0.95\t0.90\t0.80\t250\n"
                "P001\tP001_N\tnormal\t100\t97\t0.97\t0.95\t0.12\t60\t55\t0.98\t0.94\t0.90\t0.82\t0.70\t245\n",
                encoding="utf-8",
            )
            import sys
            argv = sys.argv
            sys.argv = [
                "generate_multiqc_stub.py",
                "--sample-qc", str(sample_qc),
                "--output", str(output),
            ]
            try:
                multiqc_main()
            finally:
                sys.argv = argv
            html = output.read_text(encoding="utf-8")
            self.assertIn("Quality Control Dashboard", html)
            self.assertIn("Median Coverage", html)
            self.assertIn("P001_T", html)

    def test_execution_report_generation(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            temp = Path(tempdir)
            sample_qc = temp / "sample_qc.tsv"
            maf = temp / "somatic_mutations.maf"
            tmb = temp / "tmb.tsv"
            purity = temp / "purity_ploidy.tsv"
            signatures = temp / "signature_exposures.tsv"
            clonality = temp / "clonality_summary.tsv"
            drivers = temp / "driver_matrix.tsv"
            arm_matrix = temp / "arm_level_cnv_matrix.tsv"
            features = temp / "genomic_features.tsv"
            software = temp / "software_versions.tsv"
            containers = temp / "containers.tsv"
            reference = temp / "reference_manifest.tsv"
            pipeline = temp / "pipeline_parameters.yaml"
            output = temp / "execution_report.html"

            sample_qc.write_text(
                "patient_id\tsample\tsample_type\tmean_target_coverage\tpct_target_20x\tduplicate_rate\n"
                "P001\tP001_T\ttumor\t120\t0.97\t0.10\n",
                encoding="utf-8",
            )
            maf.write_text(
                "Hugo_Symbol\tChromosome\tStart_Position\tTumor_Seq_Allele2\nTP53\t1\t5\tT\nKRAS\t12\t25398284\tA\n",
                encoding="utf-8",
            )
            tmb.write_text(
                "patient_id\ttumor_sample\tqualifying_mutations\tcallable_mb\ttmb_mut_per_mb\nP001\tP001_T\t2\t30.0\t0.0667\n",
                encoding="utf-8",
            )
            purity.write_text(
                "patient_id\tpurity\tploidy\tdiplogr\nP001\t0.72\t2.8\t-0.18\n",
                encoding="utf-8",
            )
            signatures.write_text(
                "patient_id\tSBS7_UV\tSBS2_APOBEC\tmutation_count\ttop_signature\ttotal_assigned\nP001\t0.0\t1.0\t2\tSBS13\t2.0\n",
                encoding="utf-8",
            )
            clonality.write_text(
                "patient_id\tclonal_mutations\tsubclonal_mutations\tambiguous_mutations\tunknown_mutations\tfraction_subclonal\tfraction_ambiguous\tnumber_of_clusters\nP001\t1\t1\t0\t0\t0.5\t0.0\t2\n",
                encoding="utf-8",
            )
            drivers.write_text("patient_id\tdriver_TP53\tdriver_KRAS\nP001\t1\t0\n", encoding="utf-8")
            arm_matrix.write_text("patient_id\t1p\t1q\nP001\tloss\tgain\n", encoding="utf-8")
            features.write_text("patient_id\tTMB\tpurity\nP001\t0.0667\t0.72\n", encoding="utf-8")
            software.write_text("tool\tversion\tsource\nnextflow\t26.04.6\ttest\n", encoding="utf-8")
            containers.write_text("tool\timage\tsif\ndeepvariant\tdocker://example\tcontainers/sif/dv.sif\n", encoding="utf-8")
            reference.write_text("category\tkey\tpath\nreference\tfasta\treferences/GRCh38.fa\n", encoding="utf-8")
            pipeline.write_text(
                "run_date: 2026-08-20\nsamplesheet: /tmp/samplesheet.csv\nreference_config: /tmp/references.yaml\nresults_dir: /tmp/results\n",
                encoding="utf-8",
            )

            import sys
            argv = sys.argv
            sys.argv = [
                "generate_execution_report.py",
                "--sample-qc", str(sample_qc),
                "--maf", str(maf),
                "--tmb", str(tmb),
                "--purity", str(purity),
                "--signatures", str(signatures),
                "--clonality", str(clonality),
                "--drivers", str(drivers),
                "--arm-matrix", str(arm_matrix),
                "--features", str(features),
                "--software-versions", str(software),
                "--containers", str(containers),
                "--reference-manifest", str(reference),
                "--pipeline-parameters", str(pipeline),
                "--output", str(output),
            ]
            try:
                execution_report_main()
            finally:
                sys.argv = argv
            html = output.read_text(encoding="utf-8")
            self.assertIn("Tumor/Normal WES Cohort Report", html)
            self.assertIn("Somatic Landscape", html)
            self.assertIn("TP53", html)
            self.assertIn("FACETS-derived", html)


if __name__ == "__main__":
    unittest.main()
