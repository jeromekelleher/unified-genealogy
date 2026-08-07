"""Implementation functions called by the Snakemake pipeline."""

import gzip
import pathlib
import shutil
import subprocess

import pysam


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


def bgzip(source: str, output: str) -> None:
    """Repack a gzip file as block gzip."""
    output_path = pathlib.Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(source, "rb") as source_file:
        with pysam.BGZFile(str(output_path), "wb") as output_file:
            shutil.copyfileobj(source_file, output_file)


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


def encode(icf: str, vcz: str, worker_processes: int) -> None:
    """Encode bio2zarr's intermediate format as VCF Zarr."""
    vcz_path = pathlib.Path(vcz)
    vcz_path.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "vcf2zarr",
        "encode",
        "--no-progress",
        "--worker-processes",
        str(worker_processes),
        icf,
        vcz,
    ]
    subprocess.run(command, check=True)
