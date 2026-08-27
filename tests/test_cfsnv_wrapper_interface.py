import unittest

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class CfSnvWrapperInterfaceTests(unittest.TestCase):
    def test_stdprep_module_uses_r_wrapper(self) -> None:
        module_text = (REPO_ROOT / "modules" / "cfsnv_stdprep" / "main.nf").read_text(encoding="utf-8")
        self.assertIn("Rscript ${projectDir}/scripts/cfsnv_wrapper.R STDprep", module_text)
        self.assertNotIn("cfsnv STDprep", module_text)
        self.assertIn("--snp-database ${snp_database}", module_text)

    def test_cfdnaprep_module_uses_r_wrapper(self) -> None:
        module_text = (REPO_ROOT / "modules" / "cfsnv_cfdnaprep" / "main.nf").read_text(encoding="utf-8")
        self.assertIn("Rscript ${projectDir}/scripts/cfsnv_wrapper.R cfDNAprep", module_text)
        self.assertNotIn("cfsnv cfDNAprep", module_text)
        self.assertIn("--snp-database ${snp_database}", module_text)

    def test_detectmuts_module_uses_r_wrapper(self) -> None:
        module_text = (REPO_ROOT / "modules" / "ctdna_mutect2" / "main.nf").read_text(encoding="utf-8")
        self.assertIn("Rscript ${projectDir}/scripts/cfsnv_wrapper.R DetectMuts", module_text)
        self.assertNotIn("cfsnv DetectMuts", module_text)
        self.assertIn("--snp-database ${snp_database}", module_text)
        self.assertIn("--min-hold-support ${min_hold_support}", module_text)
        self.assertIn("--min-pass-support ${min_pass_support}", module_text)

    def test_main_passes_common_snps_and_hold_thresholds_to_cfsnv(self) -> None:
        workflow_text = (REPO_ROOT / "main.nf").read_text(encoding="utf-8")
        self.assertIn("ctdna_min_hold = readConfigValue(pipeline_config_path, 'ctdna.cfsnv.min_hold')", workflow_text)
        self.assertIn("tuple(meta, plasma_r1, plasma_r2, wbc_r1, wbc_r2, bait_bed, reference_path.toString(), reference_fasta, common_snps)", workflow_text)
        self.assertIn("tuple(meta, plasma_r1, plasma_r2, bait_bed, reference_path.toString(), reference_fasta, common_snps)", workflow_text)
        self.assertIn("reference_fasta, common_snps, ctdna_blocked_positions_vcf, ctdna_min_hold, ctdna_min_pass", workflow_text)

    def test_cfsnv_dockerfile_installs_r_package_toolchain(self) -> None:
        dockerfile = (REPO_ROOT / "containers" / "definitions" / "cfsnv.Dockerfile").read_text(encoding="utf-8")
        self.assertIn("micromamba create -y -n cfsnv", dockerfile)
        self.assertIn("boost-cpp", dockerfile)
        self.assertIn("GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2", dockerfile)
        self.assertIn("picard-2.18.4.jar", dockerfile)
        self.assertIn("cfSNV_0.99.0.tar.gz", dockerfile)


if __name__ == "__main__":
    unittest.main()
