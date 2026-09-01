import unittest

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class CfSnvWrapperInterfaceTests(unittest.TestCase):
    def test_stdprep_module_uses_r_wrapper(self) -> None:
        module_text = (REPO_ROOT / "modules" / "cfsnv_stdprep" / "main.nf").read_text(encoding="utf-8")
        self.assertIn("Rscript ${projectDir}/scripts/cfsnv_wrapper.R STDprep", module_text)
        self.assertNotIn("cfsnv STDprep", module_text)
        self.assertIn("--snp-database ${snp_database}", module_text)
        self.assertIn('export TMPDIR="\\$PWD/tmp"', module_text)
        self.assertIn('export CFSNV_R_LIB_ROOT="\\$PWD/r_libs"', module_text)
        self.assertIn('export CFSNV_JAVA="/opt/conda/bin/java"', module_text)
        self.assertIn('export CFSNV_PICARD_JAR="/usr/local/share/cfsnv-tools/picard.jar"', module_text)
        self.assertIn('export CFSNV_GATK_JAR="/usr/local/share/cfsnv-tools/GenomeAnalysisTK.jar"', module_text)

    def test_cfdnaprep_module_uses_r_wrapper(self) -> None:
        module_text = (REPO_ROOT / "modules" / "cfsnv_cfdnaprep" / "main.nf").read_text(encoding="utf-8")
        self.assertIn("Rscript ${projectDir}/scripts/cfsnv_wrapper.R cfDNAprep", module_text)
        self.assertNotIn("cfsnv cfDNAprep", module_text)
        self.assertIn("--snp-database ${snp_database}", module_text)
        self.assertIn('export TMPDIR="\\$PWD/tmp"', module_text)
        self.assertIn('export CFSNV_R_LIB_ROOT="\\$PWD/r_libs"', module_text)
        self.assertIn('export CFSNV_JAVA="/opt/conda/bin/java"', module_text)
        self.assertIn('export CFSNV_PICARD_JAR="/usr/local/share/cfsnv-tools/picard.jar"', module_text)
        self.assertIn('export CFSNV_GATK_JAR="/usr/local/share/cfsnv-tools/GenomeAnalysisTK.jar"', module_text)

    def test_detectmuts_module_uses_r_wrapper(self) -> None:
        module_text = (REPO_ROOT / "modules" / "ctdna_mutect2" / "main.nf").read_text(encoding="utf-8")
        self.assertIn("Rscript ${projectDir}/scripts/cfsnv_wrapper.R DetectMuts", module_text)
        self.assertNotIn("cfsnv DetectMuts", module_text)
        self.assertIn("--snp-database ${snp_database}", module_text)
        self.assertIn("--min-hold-support ${min_hold_support}", module_text)
        self.assertIn("--min-pass-support ${min_pass_support}", module_text)
        self.assertIn('export TMPDIR="\\$PWD/tmp"', module_text)
        self.assertIn('export CFSNV_R_LIB_ROOT="\\$PWD/r_libs"', module_text)
        self.assertIn('export CFSNV_JAVA="/opt/conda/bin/java"', module_text)
        self.assertIn('export CFSNV_PICARD_JAR="/usr/local/share/cfsnv-tools/picard.jar"', module_text)

    def test_main_passes_common_snps_and_hold_thresholds_to_cfsnv(self) -> None:
        workflow_text = (REPO_ROOT / "main.nf").read_text(encoding="utf-8")
        self.assertIn("ctdna_min_hold = readConfigValue(pipeline_config_path, 'ctdna.cfsnv.min_hold')", workflow_text)
        self.assertIn("reference_fasta_path, reference_fasta_index_path, reference_fasta_dict_path", workflow_text)
        self.assertIn("reference_bwa_index_amb, reference_bwa_index_ann, reference_bwa_classic_bwt, reference_bwa_index_pac, reference_bwa_classic_sa", workflow_text)
        self.assertIn("common_snps_path, common_snps_index_path", workflow_text)
        self.assertIn("ctdna_blocked_positions_vcf_path, ctdna_blocked_positions_index_path, ctdna_min_hold, ctdna_min_pass", workflow_text)

    def test_cfsnv_dockerfile_installs_r_package_toolchain(self) -> None:
        dockerfile = (REPO_ROOT / "containers" / "definitions" / "cfsnv.Dockerfile").read_text(encoding="utf-8")
        self.assertIn("micromamba create -y -n cfsnv", dockerfile)
        self.assertIn("openjdk=8", dockerfile)
        self.assertIn("micromamba install -y -n base -c conda-forge openjdk=8", dockerfile)
        self.assertIn("boost-cpp", dockerfile)
        self.assertIn("GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2", dockerfile)
        self.assertIn("releases/download/2.18.4/picard.jar", dockerfile)
        self.assertIn("org/broadinstitute/barclay/argparser/CommandLineProgramProperties.class", dockerfile)
        self.assertIn("cfSNV_0.99.0.tar.gz", dockerfile)
        self.assertIn("ENV JAVA_HOME=/opt/conda", dockerfile)

    def test_bootstrap_derives_cfsnv_blacklist_from_common_snps_and_targets(self) -> None:
        bootstrap_text = (REPO_ROOT / "scripts" / "bootstrap_references.sh").read_text(encoding="utf-8")
        self.assertIn("ensure_cfsnv_blocked_positions()", bootstrap_text)
        self.assertIn('log "Deriving cfSNV blocked-position blacklist from ${source_vcf} on targets ${targets_bed}"', bootstrap_text)
        self.assertIn('bcftools view \\', bootstrap_text)
        self.assertIn('-R "${targets_bed}" \\', bootstrap_text)
        self.assertIn('ensure_cfsnv_blocked_positions "${PROJECT_DIR}/${common_path}" "${PROJECT_DIR}/${targets_bed}" "${PROJECT_DIR}/${cfsnv_blocked_path}"', bootstrap_text)

    def test_wrapper_copies_cfsnv_to_writable_library_and_resolves_java(self) -> None:
        wrapper_text = (REPO_ROOT / "scripts" / "cfsnv_wrapper.R").read_text(encoding="utf-8")
        self.assertIn('ensure_writable_cfsnv_library <- function()', wrapper_text)
        self.assertIn('Sys.getenv("CFSNV_R_LIB_ROOT", "")', wrapper_text)
        self.assertIn('target_pkg <- file.path(lib_root, basename(source_pkg))', wrapper_text)
        self.assertIn('file.copy(source_pkg, lib_root, recursive = TRUE)', wrapper_text)
        self.assertIn('library("cfSNV", lib.loc = lib_root, character.only = TRUE)', wrapper_text)
        self.assertIn('configure_cfsnv_tmpdir <- function(package_root, tmpdir)', wrapper_text)
        self.assertIn('package_tmpdir <- file.path(extdata_dir, "tmp")', wrapper_text)
        self.assertIn('stop("cfSNV wrapper requires Unix-style symlink support for tmpdir configuration")', wrapper_text)
        self.assertIn('file.symlink(tmpdir, package_tmpdir)', wrapper_text)
        self.assertIn('resolved_link <- Sys.readlink(package_tmpdir)', wrapper_text)
        self.assertIn('stop(sprintf("failed to bind cfSNV extdata tmpdir %s to %s", package_tmpdir, tmpdir))', wrapper_text)
        self.assertIn('configure_cfsnv_tmpdir(cfsnv_package_root, tmpdir)', wrapper_text)
        self.assertIn('output_dir <- dirname(normalizePath(require_arg(args, "output"), mustWork = FALSE))', wrapper_text)
        self.assertIn('resolve_java_path <- function()', wrapper_text)
        self.assertIn('Sys.getenv("JAVA_HOME", "")', wrapper_text)
        self.assertIn('tool_path("java", default_paths = c("/opt/conda/bin/java", "/usr/bin/java"))', wrapper_text)
        self.assertIn('resolve_picard_path <- function()', wrapper_text)
        self.assertIn('resolve_gatk3_path <- function()', wrapper_text)
        self.assertIn('require_existing_path <- function(path, label)', wrapper_text)
        self.assertIn('stage_cfsnv_input_file <- function(source_path, target_path)', wrapper_text)
        self.assertIn('stage_cfsnv_detectmuts_inputs <- function(tmpdir, sample_id, tumor_bam, normal_bam, extended_bam, not_combined_bam)', wrapper_text)
        self.assertIn('sprintf("%s.paired-reads.bam", sample_id)', wrapper_text)
        self.assertIn('sprintf("%s.normal-blood.bam", sample_id)', wrapper_text)
        self.assertIn('configure_cfsnv_tmpdir(cfsnv_package_root, tmpdir)', wrapper_text)
        self.assertIn('stage_cfsnv_detectmuts_inputs(', wrapper_text)
        self.assertIn('default_paths = c("/usr/local/share/cfsnv-tools/picard.jar")', wrapper_text)
        self.assertIn('default_paths = c("/usr/local/share/cfsnv-tools/GenomeAnalysisTK.jar", "/opt/gatk3/GenomeAnalysisTK.jar")', wrapper_text)

    def test_detectmuts_module_exports_gatk3_path(self) -> None:
        module_text = (REPO_ROOT / "modules" / "ctdna_mutect2" / "main.nf").read_text(encoding="utf-8")
        self.assertIn('export CFSNV_GATK_JAR="/usr/local/share/cfsnv-tools/GenomeAnalysisTK.jar"', module_text)

    def test_cfsnv_dockerfile_pins_ml_stack_to_upstream_supported_versions(self) -> None:
        dockerfile = (REPO_ROOT / "containers" / "definitions" / "cfsnv.Dockerfile").read_text(encoding="utf-8")
        self.assertIn("numpy==1.19.5", dockerfile)
        self.assertIn("pandas==1.1.5", dockerfile)
        self.assertIn("scipy==1.5.4", dockerfile)
        self.assertIn("scikit-learn==0.24.1", dockerfile)


if __name__ == "__main__":
    unittest.main()
