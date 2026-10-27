import pathlib

import steps


configfile: "config.yaml"


data_dir = pathlib.Path("data")
vcf_dir = data_dir / "vcf"
exploded_vcf_dir = data_dir / "exploded_vcfs"
zarr_vcf_dir = data_dir / "zarr_vcfs"
recombination_map_dir = data_dir / "recomb_map"
tsinfer_dir = data_dir / "tsinfer"
# Outputs of one tsinfer run, named by chromosome and region.
run_dir = tsinfer_dir / "chr{chrom}_{start}_{end}"


wildcard_constraints:
    dataset="tgp|sgdp|ancestral_alleles",
    chrom=r"\d+",
    start=r"\d+",
    end=r"\d+",


rule all:
    input:
        expand(zarr_vcf_dir / "tgp_chr{chrom}.vcz", chrom=config["chrom"]),
        expand(zarr_vcf_dir / "sgdp_chr{chrom}.vcz", chrom=config["chrom"]),
        expand(
            zarr_vcf_dir / "ancestral_alleles_chr{chrom}.vcz",
            chrom=config["chrom"],
        ),
        recombination_map_dir / "plink.GRCh37.map.zip",
        expand(
            run_dir / "chr{chrom}_{start}_{end}.trees",
            chrom=config["chrom"],
            start=config["region"]["start"],
            end=config["region"]["end"],
        ),


rule download_tgp_vcf:
    output:
        vcf_dir
        / (
            "ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b."
            "20130502.genotypes.vcf.gz"
        ),
    params:
        base_url=config["tgp_base_url"],
    run:
        filename = pathlib.Path(output[0]).name
        url = f"{params.base_url}/{filename}"
        steps.download(url, output[0])


rule download_tgp_index:
    output:
        vcf_dir
        / (
            "ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b."
            "20130502.genotypes.vcf.gz.tbi"
        ),
    params:
        base_url=config["tgp_base_url"],
    run:
        filename = pathlib.Path(output[0]).name
        url = f"{params.base_url}/{filename}"
        steps.download(url, output[0])


rule download_sgdp_bcf:
    output:
        vcf_dir / "chr.sgdp.pub.{chrom}.bcf",
    params:
        base_url=config["sgdp_base_url"],
    run:
        filename = pathlib.Path(output[0]).name
        url = f"{params.base_url}/{filename}"
        steps.download(url, output[0])


rule download_sgdp_index:
    output:
        vcf_dir / "chr.sgdp.pub.{chrom}.bcf.csi",
    params:
        base_url=config["sgdp_base_url"],
    run:
        filename = pathlib.Path(output[0]).name
        url = f"{params.base_url}/{filename}"
        steps.download(url, output[0])


rule download_ancestral_alleles:
    output:
        vcf_dir / "homo_sapiens-chr{chrom}.vcf.gz",
    params:
        base_url=config["ancestral_allele_base_url"],
    run:
        filename = pathlib.Path(output[0]).name
        url = f"{params.base_url}/{filename}"
        steps.download(url, output[0])


rule prepare_ancestral_alleles:
    input:
        vcf_dir / "homo_sapiens-chr{chrom}.vcf.gz",
    output:
        vcf_dir / "ancestral_alleles_chr{chrom}.vcf.gz",
    run:
        steps.prepare_ancestral_alleles(input[0], output[0], wildcards.chrom)


rule index_ancestral_alleles:
    input:
        vcf_dir / "ancestral_alleles_chr{chrom}.vcf.gz",
    output:
        vcf_dir / "ancestral_alleles_chr{chrom}.vcf.gz.tbi",
    run:
        steps.index_vcf(input[0])


rule download_recomb_map:
    output:
        recombination_map_dir / "plink.GRCh37.map.zip",
    params:
        base_url=config["recomb_map_base_url"],
    run:
        filename = pathlib.Path(output[0]).name
        url = f"{params.base_url}/{filename}"
        steps.download(url, output[0])


rule explode_tgp:
    input:
        vcf=vcf_dir
        / (
            "ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b."
            "20130502.genotypes.vcf.gz"
        ),
        index=vcf_dir
        / (
            "ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b."
            "20130502.genotypes.vcf.gz.tbi"
        ),
    output:
        directory(exploded_vcf_dir / "tgp_chr{chrom}.icf"),
    threads:
        config["worker_processes"]
    run:
        steps.explode(input.vcf, output[0], threads)


rule explode_sgdp:
    input:
        bcf=vcf_dir / "chr.sgdp.pub.{chrom}.bcf",
        index=vcf_dir / "chr.sgdp.pub.{chrom}.bcf.csi",
    output:
        directory(exploded_vcf_dir / "sgdp_chr{chrom}.icf"),
    threads:
        config["worker_processes"]
    run:
        steps.explode(input.bcf, output[0], threads)


rule explode_ancestral_alleles:
    input:
        vcf=vcf_dir / "ancestral_alleles_chr{chrom}.vcf.gz",
        index=vcf_dir / "ancestral_alleles_chr{chrom}.vcf.gz.tbi",
    output:
        directory(exploded_vcf_dir / "ancestral_alleles_chr{chrom}.icf"),
    threads:
        config["worker_processes"]
    run:
        steps.explode(input.vcf, output[0], threads)


rule mkschema:
    input:
        exploded_vcf_dir / "{dataset}_chr{chrom}.icf",
    output:
        zarr_vcf_dir / "{dataset}_chr{chrom}.schema.json",
    params:
        drop_fields=lambda wildcards: config["drop_fields"][wildcards.dataset],
    run:
        steps.make_schema(input[0], output[0], params.drop_fields)


rule encode:
    input:
        icf=exploded_vcf_dir / "{dataset}_chr{chrom}.icf",
        schema=zarr_vcf_dir / "{dataset}_chr{chrom}.schema.json",
    output:
        directory(zarr_vcf_dir / "{dataset}_chr{chrom}.vcz"),
    threads: config["worker_processes"]
    run:
        steps.encode(input.icf, input.schema, output[0], threads)


rule tsinfer_config:
    input:
        tgp=zarr_vcf_dir / "tgp_chr{chrom}.vcz",
        sgdp=zarr_vcf_dir / "sgdp_chr{chrom}.vcz",
        ancestral_alleles=zarr_vcf_dir / "ancestral_alleles_chr{chrom}.vcz",
    output:
        run_dir / "config.toml",
    run:
        output_dir = pathlib.Path(output[0]).parent
        steps.write_tsinfer_config(
            output[0],
            chrom=wildcards.chrom,
            start=int(wildcards.start),
            end=int(wildcards.end),
            sources={"tgp": input.tgp, "sgdp": input.sgdp},
            ancestral_alleles=input.ancestral_alleles,
            ancestors=output_dir / "ancestors.vcz",
            match_output=output_dir / "match.trees",
            tsinfer_config=config["tsinfer"],
        )


rule infer_ancestors:
    input:
        run_dir / "config.toml",
    output:
        directory(run_dir / "ancestors.vcz"),
    log:
        run_dir / "logs" / "infer_ancestors.log",
    threads: config["tsinfer"]["threads"]
    run:
        steps.run_tsinfer("infer-ancestors", input[0], threads, log[0], [])


rule match:
    input:
        config=run_dir / "config.toml",
        ancestors=run_dir / "ancestors.vcz",
    output:
        run_dir / "match.trees",
    log:
        run_dir / "logs" / "match.log",
    threads: config["tsinfer"]["threads"]
    run:
        steps.run_tsinfer("match", input.config, threads, log[0], [])


rule post_process:
    input:
        config=run_dir / "config.toml",
        match=run_dir / "match.trees",
    output:
        run_dir / "chr{chrom}_{start}_{end}.trees",
    log:
        run_dir / "logs" / "post_process.log",
    threads: config["tsinfer"]["threads"]
    run:
        extra_args = ["--input", input.match, "--output", output[0]]
        steps.run_tsinfer("post-process", input.config, threads, log[0], extra_args)
