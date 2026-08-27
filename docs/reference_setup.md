# Reference Setup

Bootstrap the real-run reference bundle with:

```bash
export TARGETS_BED=/absolute/path/to/assay_targets.bed
export INTERVAL_LIST=/absolute/path/to/assay_targets.interval_list
export CANCER_GENE_CENSUS=/absolute/path/to/cancer_gene_census.tsv
export HOTSPOTS_TSV=/absolute/path/to/hotspots.tsv
export GENE_COORDINATES_TSV=/absolute/path/to/gene_coordinates.tsv
export CHROMOSOME_ARMS_BED=/absolute/path/to/chromosome_arms.bed
bash scripts/bootstrap_references.sh
```

The pinned source definitions live in [config/reference_sources.yaml](/home/arc85/Desktop/wes_workflow/config/reference_sources.yaml). The bootstrap process populates:

- Broad GRCh38 FASTA, `.fai`, and sequence dictionary
- BWA-MEM2 index sidecars
- GATK germline resource, panel of normals, and derived common biallelic SNP resource
- offline Ensembl VEP release `116` cache for `homo_sapiens` / `GRCh38`
- SigProfiler reference volume under `references/GRCh38/sigprofiler`
- user-supplied assay interval files and annotation TSV files copied into the workflow tree

After bootstrap, validate with:

```bash
python3 scripts/validate_reference.py --config config/references.yaml
```
