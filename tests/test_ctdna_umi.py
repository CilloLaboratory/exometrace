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
        self.assertIn("fgbio FastqToBam", module_text)
        self.assertIn("--read-structures ${read_structure_r1} ${read_structure_r2}", module_text)
        self.assertIn("fgbio ZipperBams", module_text)
        self.assertIn("--tags-to-reverse Consensus", module_text)
        self.assertIn("--tags-to-revcomp Consensus", module_text)
        self.assertIn("read_structure_r1", module_text)
        self.assertIn("read_structure_r2", module_text)

    def test_main_routes_trimmed_fastqs_only_to_cfsnv_path(self) -> None:
        workflow_text = (REPO_ROOT / "main.nf").read_text(encoding="utf-8")
        self.assertIn("template_trimmed = UMI_TEMPLATE_TRIM", workflow_text)
        self.assertIn("consensus_inputs = fastqc_results.map", workflow_text)
        self.assertIn("ctdna_read_structure_r1", workflow_text)
        self.assertIn("ctdna_read_structure_r2", workflow_text)
        self.assertIn("tuple(meta, r1, r2, bait_bed, reference_path.toString(), reference_fasta_path, reference_fasta_index_path, reference_fasta_dict_path", workflow_text)

    def test_bwa_alignment_stages_reference_sidecars(self) -> None:
        workflow_text = (REPO_ROOT / "main.nf").read_text(encoding="utf-8")
        align_module_text = (REPO_ROOT / "modules" / "align" / "main.nf").read_text(encoding="utf-8")
        consensus_module_text = (REPO_ROOT / "modules" / "umi_consensus" / "main.nf").read_text(encoding="utf-8")

        self.assertIn('reference_fasta_path = file(reference_fasta, checkIfExists: true)', workflow_text)
        self.assertIn('reference_bwa_index_bwt = file("${reference_fasta}.bwt.2bit.64", checkIfExists: true)', workflow_text)
        self.assertIn('tuple(meta, r1, r2, bait_bed, reference_path.toString(), reference_fasta_path, reference_bwa_index_0123, reference_bwa_index_amb, reference_bwa_index_ann, reference_bwa_index_bwt, reference_bwa_index_pac)', workflow_text)
        self.assertIn('path(ref_fasta), path(ref_fasta_0123), path(ref_fasta_amb), path(ref_fasta_ann), path(ref_fasta_bwt), path(ref_fasta_pac)', align_module_text)
        self.assertIn('path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict), path(ref_fasta_0123), path(ref_fasta_amb), path(ref_fasta_ann), path(ref_fasta_bwt), path(ref_fasta_pac)', consensus_module_text)

    def test_deepvariant_and_deepsomatic_stage_reference_fasta_index(self) -> None:
        workflow_text = (REPO_ROOT / "main.nf").read_text(encoding="utf-8")
        deepvariant_text = (REPO_ROOT / "modules" / "deepvariant" / "main.nf").read_text(encoding="utf-8")
        deepsomatic_text = (REPO_ROOT / "modules" / "deepsomatic" / "main.nf").read_text(encoding="utf-8")

        self.assertIn('reference_fasta_index_path = file("${reference_fasta}.fai", checkIfExists: true)', workflow_text)
        self.assertIn('tuple(meta, bam, bai, bait_bed, ref_cfg, reference_fasta_path, reference_fasta_index_path)', workflow_text)
        self.assertIn('tuple(meta, tumor_bam, tumor_bai, normal_bam, normal_bai, bait_bed, ref_cfg, reference_fasta_path, reference_fasta_index_path)', workflow_text)
        self.assertIn('path(ref_fasta), path(ref_fasta_fai)', deepvariant_text)
        self.assertIn('path(ref_fasta), path(ref_fasta_fai)', deepsomatic_text)

    def test_reference_assets_are_staged_for_all_container_modules(self) -> None:
        workflow_text = (REPO_ROOT / "main.nf").read_text(encoding="utf-8")
        alignment_qc_text = (REPO_ROOT / "modules" / "alignment_qc" / "main.nf").read_text(encoding="utf-8")
        mutect2_text = (REPO_ROOT / "modules" / "mutect2" / "main.nf").read_text(encoding="utf-8")
        ctdna_mutect2_text = (REPO_ROOT / "modules" / "ctdna_mutect2" / "main.nf").read_text(encoding="utf-8")
        cfsnv_stdprep_text = (REPO_ROOT / "modules" / "cfsnv_stdprep" / "main.nf").read_text(encoding="utf-8")
        cfsnv_cfdnaprep_text = (REPO_ROOT / "modules" / "cfsnv_cfdnaprep" / "main.nf").read_text(encoding="utf-8")
        cnvkit_text = (REPO_ROOT / "modules" / "cnvkit" / "main.nf").read_text(encoding="utf-8")
        purity_ploidy_text = (REPO_ROOT / "modules" / "purity_ploidy" / "main.nf").read_text(encoding="utf-8")
        vep_text = (REPO_ROOT / "modules" / "vep" / "main.nf").read_text(encoding="utf-8")
        driver_text = (REPO_ROOT / "modules" / "driver_annotation" / "main.nf").read_text(encoding="utf-8")
        arm_level_text = (REPO_ROOT / "modules" / "arm_level_cnv" / "main.nf").read_text(encoding="utf-8")

        self.assertIn('reference_fasta_dict_path = file(reference_fasta.replaceFirst(/\\.fa(sta)?$/, \'.dict\'), checkIfExists: true)', workflow_text)
        self.assertIn('germline_resource_index_path = file("${germline_resource}.tbi", checkIfExists: true)', workflow_text)
        self.assertIn('panel_of_normals_index_path = file("${panel_of_normals}.tbi", checkIfExists: true)', workflow_text)
        self.assertIn('vep_cache_path = file(vep_cache, checkIfExists: true)', workflow_text)
        self.assertIn('path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict), path(target_intervals)', alignment_qc_text)
        self.assertIn('path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict), path(germline_resource), path(germline_resource_tbi), path(common_snps), path(common_snps_tbi), path(pon), path(pon_tbi), path(target_intervals)', mutect2_text)
        self.assertIn('path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict), path(germline_resource), path(germline_resource_tbi), path(common_snps), path(common_snps_tbi), path(pon), path(pon_tbi), path(target_intervals)', ctdna_mutect2_text)
        self.assertIn('path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict), path(ref_fasta_0123), path(ref_fasta_amb), path(ref_fasta_ann), path(ref_fasta_bwt), path(ref_fasta_pac), path(snp_database), path(snp_database_tbi)', cfsnv_stdprep_text)
        self.assertIn('path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict), path(ref_fasta_0123), path(ref_fasta_amb), path(ref_fasta_ann), path(ref_fasta_bwt), path(ref_fasta_pac), path(snp_database), path(snp_database_tbi)', cfsnv_cfdnaprep_text)
        self.assertIn('path(ref_fasta), path(ref_fasta_fai)', cnvkit_text)
        self.assertIn('path(common_snps), path(common_snps_tbi)', purity_ploidy_text)
        self.assertIn('path(ref_fasta), path(ref_fasta_fai), path(vep_cache)', vep_text)
        self.assertIn('path(census), path(hotspots)', driver_text)
        self.assertIn('path(arms_bed)', arm_level_text)

    def test_gatk_modules_set_explicit_java_heap(self) -> None:
        mutect2_text = (REPO_ROOT / "modules" / "mutect2" / "main.nf").read_text(encoding="utf-8")
        ctdna_mutect2_text = (REPO_ROOT / "modules" / "ctdna_mutect2" / "main.nf").read_text(encoding="utf-8")
        filter_mutect_text = (REPO_ROOT / "modules" / "filter_mutect" / "main.nf").read_text(encoding="utf-8")
        alignment_qc_text = (REPO_ROOT / "modules" / "alignment_qc" / "main.nf").read_text(encoding="utf-8")
        markduplicates_text = (REPO_ROOT / "modules" / "markduplicates" / "main.nf").read_text(encoding="utf-8")

        self.assertIn('def gatk_heap_gb = Math.max(2, task.memory.toGiga().intValue() - 4)', mutect2_text)
        self.assertIn('gatk --java-options "-Xms1g -Xmx${gatk_heap_gb}g" GetPileupSummaries', mutect2_text)
        self.assertIn('gatk --java-options "-Xms1g -Xmx${gatk_heap_gb}g" GetPileupSummaries', ctdna_mutect2_text)
        self.assertIn('gatk --java-options "-Xms1g -Xmx${gatk_heap_gb}g" FilterMutectCalls', filter_mutect_text)
        self.assertIn('gatk --java-options "-Xms1g -Xmx${gatk_heap_gb}g" CollectAlignmentSummaryMetrics', alignment_qc_text)
        self.assertIn('gatk --java-options "-Xms1g -Xmx${gatk_heap_gb}g" MarkDuplicates', markduplicates_text)

    def test_umi_consensus_sets_explicit_java_heap(self) -> None:
        module_text = (REPO_ROOT / "modules" / "umi_consensus" / "main.nf").read_text(encoding="utf-8")
        self.assertIn('def fgbio_heap_gb = Math.max(4, task.memory.toGiga().intValue() - 4)', module_text)
        self.assertIn('export JAVA_TOOL_OPTIONS="-Xms1g -Xmx${fgbio_heap_gb}g"', module_text)
        self.assertIn('export TMPDIR="\\$PWD/tmp"', module_text)

    def test_umi_consensus_uses_ubam_pipeline(self) -> None:
        module_text = (REPO_ROOT / "modules" / "umi_consensus" / "main.nf").read_text(encoding="utf-8")
        dockerfile_text = (REPO_ROOT / "containers" / "definitions" / "umi_consensus.Dockerfile").read_text(encoding="utf-8")

        self.assertIn('samtools fastq ${meta.sample_id}.unmapped.bam', module_text)
        self.assertIn('bwa mem -t ${task.cpus} -p -K 150000000 -Y ${ref_fasta} -', module_text)
        self.assertIn('samtools sort --template-coordinate --threads ${sort_threads}', module_text)
        self.assertIn('fgbio GroupReadsByUmi', module_text)
        self.assertIn('fgbio CallMolecularConsensusReads', module_text)
        self.assertIn('samtools fastq ${meta.sample_id}.consensus.unmapped.bam', module_text)
        self.assertIn('bwa=0.7.17', dockerfile_text)


if __name__ == "__main__":
    unittest.main()
