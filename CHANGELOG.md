# Changelog

## 0.1.0 - 2026-08-19

- Established repository scaffold for the tumor/normal WES workflow.
- Added Nextflow DSL2 skeleton, config files, docs, validation scripts, and unit tests.
- Pinned initial container/tool versions for workflow planning and provenance.
- Implemented Phase 2 preprocessing modules for FASTQ validation, FastQC, alignment, duplicate marking, alignment QC, and cohort QC merging.
- Added local test fixtures and a passing `-stub-run` integration path.
- Implemented Phase 3 caller modules for DeepSomatic, Mutect2 filtering, DeepVariant, and somatic callset comparison.
- Implemented Phase 4 CNV modules for CNVkit, purity/ploidy summaries, allele-specific CNV export, and arm-level cohort CNV summaries.
- Implemented Phase 5 annotation modules for VEP, MAF conversion, and driver annotation.
