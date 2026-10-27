# Unified genealogy inference with tsinfer 1.0

This Snakemake pipeline will run ARG inference from the combined Simons Genome Diversity Project (SGDP) and 1000 Genomes Project (1KGP) datasets using the new approach in tsinfer 1.0 (https://github.com/tskit-dev/tsinfer). This will ultimately result in an update to the ARGs inferred for the Unified Genealogies paper (https://www.science.org/doi/10.1126/science.abi8264) that are available on Zenodo (https://zenodo.org/records/5512994). 

The datasets are aligned to the GRCh37 reference, as previously, because there is still no phased GRCh38 SGDP dataset available and much ancient DNA research is still done using GRCh37. The multi-source inference implementation in tsinfer 1.0 is incomplete, so further development will be needed before complete ARGs can be inferred. As such, the pipeline is a work-in-progress.

## Running

```bash
uv sync
uv run snakemake --cores 4
```

The chromosomes, the region to infer and the tsinfer settings are set in
`config.yaml`. Outputs for a region are written to
`data/tsinfer/chr{chrom}_{start}_{end}/`: the generated tsinfer `config.toml`,
the ancestors (`ancestors.vcz`), the raw match output (`match.trees`), the final
post-processed tree sequence (`chr{chrom}_{start}_{end}.trees`) and the tsinfer
logs (`logs/`).

The Ensembl ancestral allele VCF is reduced to one SNV record per position
before conversion, because tsinfer needs at most one ancestral state per site.
Large per-call fields that tsinfer doesn't read (`drop_fields` in `config.yaml`)
are left out of the VCF Zarr encoding.

