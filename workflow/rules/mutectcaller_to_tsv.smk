__author__ = "Martin Rippin"
__copyright__ = "Copyright 2021, Martin Rippin"
__email__ = "martin.rippin@scilifelab.uu.se"
__license__ = "GPL-3"


rule mutectcaller_to_tsv:
    input:
        vcf="parabricks/pbrun_mutectcaller_{analysis}/{sample}_T.vep.{bed}.vcf",
    output:
        tsv=temp("tsv_files/{sample}_mutectcaller_{analysis}.{bed}.tsv"),
    params:
        analysis="{analysis}",
        sample="{sample}",
    log:
        "tsv_files/{sample}_mutectcaller_{analysis}.{bed}.tsv.log",
    benchmark:
        repeat(
            "tsv_files/{sample}_mutectcaller_{analysis}.{bed}.tsv.benchmark.tsv",
            config.get("mutectcaller_to_tsv", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("mutectcaller_to_tsv", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("mutectcaller_to_tsv", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mutectcaller_to_tsv", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("mutectcaller_to_tsv", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("mutectcaller_to_tsv", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("mutectcaller_to_tsv", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("mutectcaller_to_tsv", {}).get("container", config["default_container"])
    message:
        "{rule}: select relevant info from {input.vcf} and compile into tsv"
    script:
        "../scripts/mutectcaller_to_tsv.py"
