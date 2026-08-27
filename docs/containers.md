# Container Strategy

Container paths and pinned versions are declared in `config/containers.yaml`.

Pins were selected from upstream official release pages checked on August 19, 2026. They should remain fixed unless deliberately updated with corresponding changelog and validation work.

Planned build location:

- definitions: `containers/definitions/`
- built SIF images: `containers/sif/`

ctDNA-specific container build assets now live in:

- `containers/definitions/cfsnv.Dockerfile`
- `containers/definitions/umi_consensus.Dockerfile`
- `scripts/build_ctdna_containers.sh`
