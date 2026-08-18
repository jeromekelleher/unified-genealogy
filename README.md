# Unified genealogy inference with tsinfer 1.0

This Snakemake pipeline will run ARG inference from the combined Simons Genome Diversity Project (SGDP) and 1000 Genomes Project (1KGP) datasets using the new approach in tsinfer 1.0 (https://github.com/tskit-dev/tsinfer). This will ultimately result in an update to the ARGs inferred for the Unified Genealogies paper (https://www.science.org/doi/10.1126/science.abi8264) that are available on Zenodo (https://zenodo.org/records/5512994). 

The datasets are aligned to the GRCh37 reference, as previously, because there is still no phased GRCh38 SGDP dataset available and much ancient DNA research is still done using GRCh37. The multi-source inference implementation in tsinfer 1.0 is incomplete, so further development will be needed before complete ARGs can be inferred. As such, the pipeline is a work-in-progress.
