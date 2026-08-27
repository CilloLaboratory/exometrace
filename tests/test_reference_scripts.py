import tempfile
import unittest

from pathlib import Path

from scripts.validate_reference import main as validate_reference_main


def write_file(path: Path, content: str = "x") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


class ReferenceValidationTests(unittest.TestCase):
    def build_config(self, root: Path) -> Path:
        config = root / "config" / "references.yaml"
        config.parent.mkdir(parents=True, exist_ok=True)
        config.write_text(
            "\n".join(
                [
                    "reference:",
                    "  build: GRCh38",
                    "  contig_style: chr",
                    "  root: references/GRCh38",
                    "  fasta: references/GRCh38/fasta/GRCh38.fa",
                    "  fasta_index: references/GRCh38/fasta/GRCh38.fa.fai",
                    "  sequence_dictionary: references/GRCh38/fasta/GRCh38.dict",
                    "  bwa_index_prefix: references/GRCh38/fasta/GRCh38.fa",
                    "",
                    "intervals:",
                    "  targets_bed: references/GRCh38/intervals/exome_targets.bed",
                    "  interval_list: references/GRCh38/intervals/exome_targets.interval_list",
                    "  callable_regions: references/GRCh38/intervals/callable_regions.bed",
                    "",
                    "gatk:",
                    "  germline_resource: references/GRCh38/gatk/af-only-gnomad.vcf.gz",
                    "  common_snps: references/GRCh38/gatk/common_biallelic_snps.vcf.gz",
                    "  panel_of_normals: references/GRCh38/gatk/panel_of_normals.vcf.gz",
                    "",
                    "vep:",
                    "  cache_dir: references/GRCh38/vep/cache",
                    "",
                    "sigprofiler:",
                    "  volume_dir: references/GRCh38/sigprofiler",
                    "",
                    "annotations:",
                    "  cancer_gene_census: references/GRCh38/annotations/cancer_gene_census.tsv",
                    "  hotspots: references/GRCh38/annotations/hotspots.tsv",
                    "  gene_coordinates: references/GRCh38/annotations/gene_coordinates.tsv",
                    "  chromosome_arms: references/GRCh38/annotations/chromosome_arms.bed",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        return config

    def populate_valid_tree(self, root: Path) -> None:
        ref_root = root / "references" / "GRCh38"
        fasta = ref_root / "fasta" / "GRCh38.fa"
        write_file(fasta, ">chr1\nACGT\n")
        write_file(ref_root / "fasta" / "GRCh38.fa.fai")
        write_file(ref_root / "fasta" / "GRCh38.dict")
        for suffix in [".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"]:
            write_file(Path(f"{fasta}{suffix}"))
        write_file(ref_root / "intervals/exome_targets.bed", "\n".join(f"chr1\t{i}\t{i+1}" for i in range(20)) + "\n")
        write_file(
            ref_root / "intervals/exome_targets.interval_list",
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100\n" + "\n".join(f"chr1\t{i}\t{i+1}\t+\ttarget_{i}" for i in range(20)) + "\n",
        )
        write_file(ref_root / "intervals/callable_regions.bed", "\n".join(f"chr1\t{i}\t{i+1}" for i in range(20)) + "\n")
        write_file(ref_root / "annotations/cancer_gene_census.tsv", "gene\tsource\n" + "\n".join(f"GENE{i}\tCGC" for i in range(20)) + "\n")
        write_file(ref_root / "annotations/hotspots.tsv", "gene\tprotein_change\n" + "\n".join(f"GENE{i}\tp.V{i}A" for i in range(20)) + "\n")
        write_file(ref_root / "annotations/gene_coordinates.tsv", "gene\tchrom\tstart\tend\n" + "\n".join(f"GENE{i}\tchr1\t{i}\t{i+10}" for i in range(20)) + "\n")
        write_file(ref_root / "annotations/chromosome_arms.bed", "\n".join(f"chr1\t{i}\t{i+1}\t1p" for i in range(20)) + "\n")
        for relative in [
            "gatk/af-only-gnomad.vcf.gz",
            "gatk/common_biallelic_snps.vcf.gz",
            "gatk/panel_of_normals.vcf.gz",
        ]:
            write_file(ref_root / relative)
            write_file(ref_root / f"{relative}.tbi")
        write_file(ref_root / "vep/cache/homo_sapiens_vep_116_GRCh38/README", "cache\n")
        write_file(ref_root / "sigprofiler/tsb/GRCh38/README", "sig\n")

    def run_validator(self, config: Path) -> int:
        import sys

        argv = sys.argv
        sys.argv = ["validate_reference.py", "--config", str(config)]
        try:
            return validate_reference_main()
        finally:
            sys.argv = argv

    def run_validator_with_pipeline(self, config: Path, pipeline_config: Path) -> int:
        import sys

        argv = sys.argv
        sys.argv = ["validate_reference.py", "--config", str(config), "--pipeline-config", str(pipeline_config)]
        try:
            return validate_reference_main()
        finally:
            sys.argv = argv

    def test_validator_accepts_real_run_tree(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            config = self.build_config(root)
            self.populate_valid_tree(root)
            self.assertEqual(self.run_validator(config), 0)

    def test_validator_requires_bwa_sidecars(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            config = self.build_config(root)
            self.populate_valid_tree(root)
            (root / "references" / "GRCh38" / "fasta" / "GRCh38.fa.pac").unlink()
            self.assertEqual(self.run_validator(config), 1)

    def test_validator_requires_populated_vep_cache(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            config = self.build_config(root)
            self.populate_valid_tree(root)
            cache_dir = root / "references" / "GRCh38" / "vep" / "cache"
            for child in cache_dir.iterdir():
                if child.is_file():
                    child.unlink()
                else:
                    for nested in child.rglob("*"):
                        if nested.is_file():
                            nested.unlink()
                    child.rmdir()
            self.assertEqual(self.run_validator(config), 1)

    def test_validator_requires_sigprofiler_assets(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            config = self.build_config(root)
            self.populate_valid_tree(root)
            sig_root = root / "references" / "GRCh38" / "sigprofiler"
            for nested in sorted(sig_root.rglob("*"), reverse=True):
                if nested.is_file():
                    nested.unlink()
                elif nested.is_dir():
                    nested.rmdir()
            sig_root.mkdir(parents=True, exist_ok=True)
            self.assertEqual(self.run_validator(config), 1)

    def test_validator_rejects_fixture_grade_annotation_tables(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            config = self.build_config(root)
            self.populate_valid_tree(root)
            write_file(root / "references" / "GRCh38" / "annotations" / "hotspots.tsv", "gene\tprotein_change\nBRAF\tp.V600E\n")
            self.assertEqual(self.run_validator(config), 1)

    def test_ctdna_mode_requires_blocked_positions_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            config = self.build_config(root)
            self.populate_valid_tree(root)
            config.write_text(
                config.read_text(encoding="utf-8")
                + "cfsnv:\n"
                + "  blocked_positions_vcf: references/GRCh38/cfsnv/blocked_positions.vcf.gz\n"
                + "  blocked_positions_index: references/GRCh38/cfsnv/blocked_positions.vcf.gz.tbi\n",
                encoding="utf-8",
            )
            pipeline_config = root / "config" / "default.yaml"
            pipeline_config.write_text("analysis:\n  mode: ctdna_umi\n", encoding="utf-8")
            self.assertEqual(self.run_validator_with_pipeline(config, pipeline_config), 1)
