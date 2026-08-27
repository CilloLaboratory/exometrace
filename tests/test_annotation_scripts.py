import unittest

from pathlib import Path
import csv
import tempfile

from scripts.annotate_drivers import main as driver_main
from scripts.ctdna_vcf_to_maf import main as ctdna_maf_main
from scripts.compare_ctdna_callsets import main as ctdna_compare_main
from scripts.vep_to_maf import main as maf_main


class AnnotationFixturesTests(unittest.TestCase):
    def test_vep_to_maf_conversion(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            output = Path(tempdir) / "out.maf"
            import sys
            argv = sys.argv
            sys.argv = [
                "vep_to_maf.py",
                "--input", "tests/test_data/annotation/deepsomatic_annotated.vcf.gz",
                "--comparison", "tests/test_data/annotation/comparison.tsv",
                "--output", str(output),
            ]
            try:
                maf_main()
            finally:
                sys.argv = argv
            with output.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["Hugo_Symbol"], "TP53")
            self.assertEqual(rows[0]["mutect2_support"], "1")

    def test_driver_annotation(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            maf_path = Path(tempdir) / "input.maf"
            maf_path.write_text(
                "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele2\tTumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\tVariant_Classification\tVariant_Type\tHGVSc\tHGVSp\tt_depth\tt_ref_count\tt_alt_count\tn_depth\tn_ref_count\tn_alt_count\ttumor_vaf\tcaller\tmutect2_support\tcaller_count\n"
                "TP53\t1\t5\t5\tA\tT\tP001_T\tP001_N\tMissense_Mutation\tSNP\tc.215C>T\tp.Pro72Leu\t30\t20\t10\t28\t28\t0\t0.333333\tdeepsomatic\t1\t2\n",
                encoding="utf-8",
            )
            long_path = Path(tempdir) / "drivers_long.tsv"
            matrix_path = Path(tempdir) / "driver_matrix.tsv"
            import sys
            argv = sys.argv
            sys.argv = [
                "annotate_drivers.py",
                "--maf", str(maf_path),
                "--census", "tests/test_data/refs/GRCh38/annotations/cancer_gene_census.tsv",
                "--hotspots", "tests/test_data/refs/GRCh38/annotations/hotspots.tsv",
                "--output-long", str(long_path),
                "--output-matrix", str(matrix_path),
            ]
            try:
                driver_main()
            finally:
                sys.argv = argv
            self.assertTrue(long_path.exists())
            self.assertTrue(matrix_path.exists())
            with long_path.open(encoding="utf-8") as handle:
                long_rows = list(csv.DictReader(handle, delimiter="\t"))
            with matrix_path.open(encoding="utf-8") as handle:
                matrix_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(long_rows[0]["patient_id"], "P001_T")
            self.assertEqual(matrix_rows[0]["patient_id"], "P001_T")

    def test_ctdna_compare_labels_support_classes(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            temp_root = Path(tempdir)
            cfsnv_vcf = temp_root / "cfsnv.vcf"
            cfsnv_vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPLASMA\tWBC\n"
                "chr1\t5\t.\tA\tT\t60\tPASS\tUMI_FAMILY_COUNT=4;UMI_MAX_FAMILY_SIZE=3\tGT:DP:AD\t0/1:28:20,8\t0/0:32:32,0\n"
                "chr1\t8\t.\tC\tG\t60\tPASS\tUMI_FAMILY_COUNT=2;UMI_MAX_FAMILY_SIZE=2\tGT:DP:AD\t0/1:20:18,2\t0/0:30:30,0\n",
                encoding="utf-8",
            )
            mutect_vcf = temp_root / "mutect2.vcf"
            mutect_vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPLASMA\tWBC\n"
                "chr1\t5\t.\tA\tT\t60\tPASS\t.\tGT:DP:AD\t0/1:31:20,11\t0/0:29:29,0\n"
                "chr1\t11\t.\tG\tC\t45\tPASS\t.\tGT:DP:AD\t0/1:18:15,3\t0/0:20:20,0\n",
                encoding="utf-8",
            )
            output = temp_root / "compare.tsv"
            import sys
            argv = sys.argv
            sys.argv = ["compare_ctdna_callsets.py", "--cfsnv", str(cfsnv_vcf), "--mutect2", str(mutect_vcf), "--output", str(output)]
            try:
                ctdna_compare_main()
            finally:
                sys.argv = argv
            with output.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            by_key = {(row["CHROM"], row["POS"], row["REF"], row["ALT"]): row for row in rows}
            self.assertEqual(by_key[("chr1", "5", "A", "T")]["support_class"], "shared")
            self.assertEqual(by_key[("chr1", "8", "C", "G")]["support_class"], "cfSNV-only")
            self.assertEqual(by_key[("chr1", "11", "G", "C")]["support_class"], "Mutect2-only")

    def test_ctdna_maf_emits_both_call_tiers_and_fields(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            temp_root = Path(tempdir)
            annotated_vcf = temp_root / "cfsnv.vep.vcf"
            annotated_vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|HGVSc|HGVSp\">\n"
                "##INFO=<ID=UMI_FAMILY_COUNT,Number=1,Type=Integer,Description=\"Consensus family support\">\n"
                "##INFO=<ID=UMI_MAX_FAMILY_SIZE,Number=1,Type=Integer,Description=\"Largest observed family size\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPLASMA\tWBC\n"
                "chr1\t5\t.\tA\tT\t60\tPASS\tUMI_FAMILY_COUNT=4;UMI_MAX_FAMILY_SIZE=3;CSQ=T|missense_variant|MODERATE|TP53|ENSG|Transcript|ENST|protein_coding|ENST:c.215C>T|ENSP:p.Pro72Leu\tGT:DP:AD\t0/1:28:20,8\t0/0:32:32,0\n"
                "chr1\t8\t.\tC\tG\t60\tPASS\tUMI_FAMILY_COUNT=1;UMI_MAX_FAMILY_SIZE=1;CSQ=G|missense_variant|MODERATE|EGFR|ENSG|Transcript|ENST|protein_coding|ENST:c.2573T>G|ENSP:p.Leu858Arg\tGT:DP:AD\t0/1:20:18,2\t0/0:30:30,0\n",
                encoding="utf-8",
            )
            comparison = temp_root / "compare.tsv"
            comparison.write_text(
                "CHROM\tPOS\tREF\tALT\tcfsnv\tmutect2\tsupport_class\tconsensus_t_depth\tconsensus_alt_count\tumi_family_count\tumi_max_family_size\twbc_alt_count\n"
                "chr1\t5\tA\tT\t1\t1\tshared\t31\t11\t4\t3\t0\n"
                "chr1\t8\tC\tG\t1\t0\tcfSNV-only\t0\t0\t1\t1\t0\n",
                encoding="utf-8",
            )
            umi_qc = temp_root / "umi_qc.tsv"
            umi_qc.write_text(
                "sample_id\tumi_family_count\tumi_max_family_size\tconsensus_reads\tstatus\n"
                "PLASMA\t10\t4\t100\tPASS\n",
                encoding="utf-8",
            )
            hs_output = temp_root / "hs.maf.tsv"
            hc_output = temp_root / "hc.maf.tsv"
            import sys
            argv = sys.argv
            sys.argv = [
                "ctdna_vcf_to_maf.py",
                "--input", str(annotated_vcf),
                "--comparison", str(comparison),
                "--umi-qc", str(umi_qc),
                "--high-sensitivity-output", str(hs_output),
                "--high-confidence-output", str(hc_output),
                "--min-family-support", "3",
            ]
            try:
                ctdna_maf_main()
            finally:
                sys.argv = argv
            with hs_output.open(encoding="utf-8") as handle:
                hs_rows = list(csv.DictReader(handle, delimiter="\t"))
            with hc_output.open(encoding="utf-8") as handle:
                hc_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(hs_rows), 2)
            self.assertEqual(len(hc_rows), 1)
            self.assertEqual(hc_rows[0]["call_tier"], "high_confidence")
            self.assertEqual(hs_rows[0]["primary_caller"], "cfSNV")
            self.assertIn("consensus_t_depth", hs_rows[0])
            self.assertIn("umi_family_count", hs_rows[0])
            self.assertEqual(hs_rows[0]["mutect2_support"], "1")


if __name__ == "__main__":
    unittest.main()
