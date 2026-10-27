"""Implementation functions called by the Snakemake pipeline."""

import json
import pathlib
import subprocess

import pysam
import pysam.bcftools
import tomli_w


def download(url: str, output: str) -> None:
    """Download ``url`` to the path managed by Snakemake."""
    output_path = pathlib.Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    partial_path = output_path.with_suffix(f"{output_path.suffix}.part")
    command = [
        "curl",
        "--fail",
        "--location",
        "--continue-at",
        "-",
        "--output",
        partial_path,
        url,
    ]
    subprocess.run(command, check=True)
    partial_path.replace(output_path)


def prepare_ancestral_alleles(source: str, output: str, chrom: str) -> None:
    """Reduce the Ensembl variation VCF to one SNV record per position.

    tsinfer requires at most one ancestral-state record at each site. The
    Ensembl VCF often has several SNV records at a position (one per rsID,
    some with no ALT allele), which all agree on REF and AA. We keep SNVs that
    have an AA annotation and merge records at the same position into one
    multi-allelic record.

    The Ensembl header declares no contigs, which bcftools needs to write BCF,
    so we first add a ``##contig`` line for ``chrom``.
    """
    output_path = pathlib.Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    header_path = output_path.with_name(f"{output_path.name}.header.txt")
    with_contig_path = output_path.with_name(f"{output_path.name}.contig.bcf")
    snvs_path = output_path.with_name(f"{output_path.name}.snvs.bcf")
    header_path.write_text(f"##contig=<ID={chrom}>\n")
    pysam.bcftools.annotate(
        "--header-lines",
        str(header_path),
        "--output-type",
        "u",
        "--output",
        str(with_contig_path),
        source,
        catch_stdout=False,
    )
    pysam.bcftools.view(
        "--types",
        "snps",
        "--include",
        'INFO/AA!="."',
        "--output-type",
        "u",
        "--output",
        str(snvs_path),
        str(with_contig_path),
        catch_stdout=False,
    )
    pysam.bcftools.norm(
        "--multiallelics",
        "+snps",
        "--output-type",
        "z",
        "--output",
        str(output_path),
        str(snvs_path),
        catch_stdout=False,
    )
    for path in [header_path, with_contig_path, snvs_path]:
        path.unlink()


def index_vcf(vcf: str) -> None:
    """Create a tabix index alongside a VCF."""
    pysam.tabix_index(vcf, preset="vcf")


def explode(vcf: str, icf: str, worker_processes: int) -> None:
    """Explode a VCF or BCF into bio2zarr's intermediate format."""
    icf_path = pathlib.Path(icf)
    icf_path.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "vcf2zarr",
        "explode",
        "--no-progress",
        "--worker-processes",
        str(worker_processes),
        vcf,
        icf,
    ]
    subprocess.run(command, check=True)


def make_schema(icf: str, schema: str, drop_fields: list[str]) -> None:
    """Write a VCF Zarr encoding schema that omits ``drop_fields``.

    Large per-call fields that tsinfer never reads (e.g. genotype likelihoods)
    dominate encode time and output size, so we leave them out.
    """
    result = subprocess.run(
        ["vcf2zarr", "mkschema", icf], check=True, capture_output=True, text=True
    )
    schema_dict = json.loads(result.stdout)
    field_names = {field["name"] for field in schema_dict["fields"]}
    unknown = set(drop_fields) - field_names
    if len(unknown) > 0:
        raise ValueError(f"Fields to drop not in schema: {sorted(unknown)}")
    kept_fields = [f for f in schema_dict["fields"] if f["name"] not in drop_fields]
    schema_dict["fields"] = kept_fields
    pathlib.Path(schema).write_text(json.dumps(schema_dict, indent=2))


def encode(icf: str, schema: str, vcz: str, worker_processes: int) -> None:
    """Encode bio2zarr's intermediate format as VCF Zarr using ``schema``."""
    vcz_path = pathlib.Path(vcz)
    vcz_path.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "vcf2zarr",
        "encode",
        "--no-progress",
        "--worker-processes",
        str(worker_processes),
        "--schema",
        schema,
        icf,
        vcz,
    ]
    subprocess.run(command, check=True)


def write_tsinfer_config(
    output: str,
    chrom: str,
    start: int,
    end: int,
    sources: dict[str, str],
    ancestral_alleles: str,
    ancestors: pathlib.Path,
    match_output: pathlib.Path,
    tsinfer_config: dict,
) -> None:
    """Write the tsinfer TOML config for inferring ``chrom:start-end``.

    ``sources`` maps each source name to its VCZ path. Every source is a sample
    source in the output tree sequence and contributes to ancestor inference.
    Paths are written as absolute paths because tsinfer resolves relative
    paths against the config file's directory.
    """
    region = f"{chrom}:{start}-{end}"
    site_include = tsinfer_config["site_include"]
    source_tables = [
        {
            "name": name,
            "path": str(pathlib.Path(path).resolve()),
            "regions": region,
            "include": site_include,
        }
        for name, path in sources.items()
    ]
    match_sources = {
        "ancestors": {"node_flags": 0, "create_individuals": False},
    }
    for name in sources:
        match_sources[name] = {"node_flags": 1, "create_individuals": True}
    toml_dict = {
        "source": source_tables,
        "ancestral_state": {
            "path": str(pathlib.Path(ancestral_alleles).resolve()),
            "field": "variant_AA",
            "regions": region,
            "include": site_include,
        },
        "ancestors": [
            {
                "name": "ancestors",
                "path": str(pathlib.Path(ancestors).resolve()),
                "sources": list(sources),
                "max_gap_length": tsinfer_config["max_gap_length"],
                "samples_chunk_size": tsinfer_config["samples_chunk_size"],
                "variants_chunk_size": tsinfer_config["variants_chunk_size"],
            }
        ],
        "match": {
            "output": str(pathlib.Path(match_output).resolve()),
            "sources": match_sources,
        },
        "post_process": {"split_ultimate": True, "erase_flanks": True},
    }
    output_path = pathlib.Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as f:
        tomli_w.dump(toml_dict, f)


def run_tsinfer(
    subcommand: str,
    config_path: str,
    threads: int,
    log_file: str,
    extra_args: list[str],
) -> None:
    """Run a ``tsinfer`` CLI subcommand on ``config_path``."""
    command = [
        "tsinfer",
        subcommand,
        config_path,
        "--threads",
        str(threads),
        "--verbose",
        "--log-file",
        log_file,
        *extra_args,
    ]
    subprocess.run(command, check=True)
