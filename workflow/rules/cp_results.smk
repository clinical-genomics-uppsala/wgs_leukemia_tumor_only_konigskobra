__author__ = "Arielle R Munters, Nina Hollfelder"
__copyright__ = "Copyright 2022, Arielle R Munters"
__email__ = "arielle.munters@scilifelab.uu.se, nina.hollfelder@scilifelab.uu.se"
__license__ = "GPL-3"


localrules:
    cp_cnvkit_t,
    cp_tbi_t,
    cp_all_vcf,
    cp_all_tbi,
    cp_aml_vcf,
    cp_aml_tbi,
    cp_cram,
    cp_crai,
    cp_tsv_mutectcaller_all,
    cp_tsv_manta_all,
    cp_tsv_manta_no_diagnosis,
    cp_tsv_mutectcaller_aml,
    cp_tsv_manta_aml,
    cp_manta_vcf,
    cp_manta_tbi,
    cp_manta_all_vcf,
    cp_manta_all_tbi,
    cp_manta_aml_vcf,
    cp_manta_aml_tbi,
    cp_pindel_vcf,
    cp_pindel_tbi,
    cp_cnvkit_vcf,
    cp_cnvkit_tbi,
    cp_cnvkit_table,
    cp_cnvkit_diagram,
    cp_cnvkit_scatter_whole,
    cp_cnvkit_scatter,
    cp_multiqc,
    cp_spring_archive,


rule cp_vcf_t:
    input:
        "parabricks/pbrun_mutectcaller_t/{sample}_T.vep.vcf.gz",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}_T.vep.vcf.gz",
    threads: config.get("cp_vcf_t", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_vcf_t", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_vcf_t", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_vcf_t", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_vcf_t", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_vcf_t", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_vcf_t", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_tbi_t:
    input:
        "parabricks/pbrun_mutectcaller_t/{sample}_T.vep.vcf.gz.tbi",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}_T.vep.vcf.gz.tbi",
    threads: config.get("cp_tbi_t", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_tbi_t", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_tbi_t", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_tbi_t", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_tbi_t", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_tbi_t", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_tbi_t", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_all_vcf:
    input:
        "parabricks/pbrun_mutectcaller_t/{sample}_T.vep.all.vcf.gz",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}_T.vep.all.vcf.gz",
    threads: config.get("cp_all_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_all_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_all_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_all_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_all_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_all_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_all_vcf", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_all_tbi:
    input:
        "parabricks/pbrun_mutectcaller_t/{sample}_T.vep.all.vcf.gz.tbi",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}_T.vep.all.vcf.gz.tbi",
    threads: config.get("cp_all_tbi", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_all_tbi", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_all_tbi", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_all_tbi", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_all_tbi", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_all_tbi", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_all_tbi", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_aml_vcf:
    input:
        "parabricks/pbrun_mutectcaller_t/{sample}_T.vep.aml.vcf.gz",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}_T.vep.aml.vcf.gz",
    threads: config.get("cp_aml_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_aml_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_aml_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_aml_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_aml_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_aml_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_aml_vcf", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_aml_tbi:
    input:
        "parabricks/pbrun_mutectcaller_t/{sample}_T.vep.aml.vcf.gz.tbi",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}_T.vep.aml.vcf.gz.tbi",
    threads: config.get("cp_aml_tbi", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_aml_tbi", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_aml_tbi", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_aml_tbi", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_aml_tbi", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_aml_tbi", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_aml_tbi", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_cram:
    input:
        "compression/crumble/{sample}_T.crumble.cram",
    output:
        "Results/{project}/{sample}/Cram/{sample}_T.crumble.cram",
    threads: config.get("cp_cram", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_cram", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_cram", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_cram", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_cram", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_cram", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_cram", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_crai:
    input:
        "compression/crumble/{sample}_T.crumble.cram.crai",
    output:
        "Results/{project}/{sample}/Cram/{sample}_T.crumble.cram.crai",
    threads: config.get("cp_crai", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_crai", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_crai", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_crai", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_crai", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_crai", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_crai", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_tsv_mutectcaller_all:
    input:
        "tsv_files/{sample}_mutectcaller_t.all.tsv",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}_mutectcaller_T.all.tsv",
    threads: config.get("cp_tsv_mutectcaller_all", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_tsv_mutectcaller_all", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_tsv_mutectcaller_all", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_tsv_mutectcaller_all", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_tsv_mutectcaller_all", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_tsv_mutectcaller_all", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_tsv_mutectcaller_all", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_tsv_manta_all:
    input:
        "tsv_files/{sample}_manta_t.all.tsv",
    output:
        "Results/{project}/{sample}/SV/{sample}_manta_T.all.tsv",
    threads: config.get("cp_tsv_manta_all", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_tsv_manta_all", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_tsv_manta_all", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_tsv_manta_all", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_tsv_manta_all", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_tsv_manta_all", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_tsv_manta_all", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_tsv_manta_no_diagnosis:
    input:
        "tsv_files/{sample}_manta_t.tsv",
    output:
        "Results/{project}/{sample}/SV/{sample}_manta_T.tsv",
    threads: config.get("cp_tsv_manta_no_diagnosis", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_tsv_manta_no_diagnosis", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_tsv_manta_no_diagnosis", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_tsv_manta_no_diagnosis", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_tsv_manta_no_diagnosis", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_tsv_manta_no_diagnosis", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_tsv_manta_no_diagnosis", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_tsv_mutectcaller_aml:
    input:
        "tsv_files/{sample}_mutectcaller_t.aml.tsv",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}_mutectcaller_T.aml.tsv",
    threads: config.get("cp_tsv_mutectcaller_aml", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_tsv_mutectcaller_aml", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_tsv_mutectcaller_aml", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_tsv_mutectcaller_aml", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_tsv_mutectcaller_aml", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_tsv_mutectcaller_aml", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_tsv_mutectcaller_aml", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_tsv_manta_aml:
    input:
        "tsv_files/{sample}_manta_t.aml.tsv",
    output:
        "Results/{project}/{sample}/SV/{sample}_manta_T.aml.tsv",
    threads: config.get("cp_tsv_manta_aml", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_tsv_manta_aml", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_tsv_manta_aml", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_tsv_manta_aml", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_tsv_manta_aml", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_tsv_manta_aml", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_tsv_manta_aml", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_manta_vcf:
    input:
        "cnv_sv/manta_run_workflow_t/{sample}.ssa.vcf.gz",
    output:
        "Results/{project}/{sample}/SV/{sample}_T.ssa.vcf.gz",
    threads: config.get("cp_manta_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_manta_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_manta_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_manta_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_manta_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_manta_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_manta_vcf", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_manta_tbi:
    input:
        "cnv_sv/manta_run_workflow_t/{sample}.ssa.vcf.gz.tbi",
    output:
        "Results/{project}/{sample}/SV/{sample}_T.ssa.vcf.gz.tbi",
    threads: config.get("cp_manta_tbi", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_manta_tbi", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_manta_tbi", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_manta_tbi", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_manta_tbi", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_manta_tbi", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_manta_tbi", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_manta_all_vcf:
    input:
        "cnv_sv/manta_run_workflow_t/{sample}.ssa.all.vcf.gz",
    output:
        "Results/{project}/{sample}/SV/{sample}_T.ssa.all.vcf.gz",
    threads: config.get("cp_manta_all_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_manta_all_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_manta_all_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_manta_all_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_manta_all_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_manta_all_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_manta_all_vcf", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"

        
rule cp_manta_all_tbi:
    input:
        "cnv_sv/manta_run_workflow_t/{sample}.ssa.all.vcf.gz.tbi",
    output:
        "Results/{project}/{sample}/SV/{sample}_T.ssa.all.vcf.gz.tbi",
    threads: config.get("cp_manta_all_tbi", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_manta_all_tbi", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_manta_all_tbi", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_manta_all_tbi", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_manta_all_tbi", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_manta_all_tbi", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_manta_all_tbi", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_manta_aml_vcf:
    input:
        "cnv_sv/manta_run_workflow_t/{sample}.ssa.aml.vcf.gz",
    output:
        "Results/{project}/{sample}/SV/{sample}_T.ssa.aml.vcf.gz",
    threads: config.get("cp_manta_aml_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_manta_aml_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_manta_aml_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_manta_aml_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_manta_aml_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_manta_aml_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_manta_aml_vcf", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_manta_aml_tbi:
    input:
        "cnv_sv/manta_run_workflow_t/{sample}.ssa.aml.vcf.gz.tbi",
    output:
        "Results/{project}/{sample}/SV/{sample}_T.ssa.aml.vcf.gz.tbi",
    threads: config.get("cp_manta_aml_tbi", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_manta_aml_tbi", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_manta_aml_tbi", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_manta_aml_tbi", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_manta_aml_tbi", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_manta_aml_tbi", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_manta_aml_tbi", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_pindel_vcf:
    input:
        "cnv_sv/pindel/{sample}.vcf.gz",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}.pindel.vcf.gz",
    threads: config.get("cp_pindel_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_pindel_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_pindel_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_pindel_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_pindel_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_pindel_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_pindel_vcf", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_pindel_tbi:
    input:
        "cnv_sv/pindel/{sample}.vcf.gz.tbi",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}.pindel.vcf.gz.tbi",
    threads: config.get("cp_pindel_tbi", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_pindel_tbi", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_pindel_tbi", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_pindel_tbi", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_pindel_tbi", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_pindel_tbi", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_pindel_tbi", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_cnvkit_vcf:
    input:
        "cnv_sv/cnvkit_vcf/{sample}_T.vcf.gz",
    output:
        "Results/{project}/{sample}/CNV/{sample}_T.vcf.gz",
    threads: config.get("cp_cnvkit_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_cnvkit_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_cnvkit_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_cnvkit_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_cnvkit_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_cnvkit_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_cnvkit_vcf", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_cnvkit_tbi:
    input:
        "cnv_sv/cnvkit_vcf/{sample}_T.vcf.gz.tbi",
    output:
        "Results/{project}/{sample}/CNV/{sample}_T.vcf.gz.tbi",
    threads: config.get("cp_cnvkit_tbi", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_cnvkit_tbi", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_cnvkit_tbi", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_cnvkit_tbi", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_cnvkit_tbi", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_cnvkit_tbi", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_cnvkit_tbi", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_cnvkit_table:
    input:
        "cnv_sv/cnvkit_table/{sample}_T.CNV.xlsx",
    output:
        "Results/{project}/{sample}/CNV/{sample}_T.CNV.xlsx",
    threads: config.get("cp_cnvkit_table", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_cnvkit_table", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_cnvkit_table", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_cnvkit_table", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_cnvkit_table", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_cnvkit_table", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_cnvkit_table", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_cnvkit_diagram:
    input:
        "cnv_sv/cnvkit_diagram/{sample}_T.png",
    output:
        "Results/{project}/{sample}/CNV/{sample}_T_diagram.png",
    threads: config.get("cp_cnvkit_diagram", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_cnvkit_diagram", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_cnvkit_diagram", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_cnvkit_diagram", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_cnvkit_diagram", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_cnvkit_diagram", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_cnvkit_diagram", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_cnvkit_scatter_whole:
    input:
        "cnv_sv/cnvkit_scatter/{sample}_T.png",
    output:
        "Results/{project}/{sample}/CNV/{sample}_T.png",
    threads: config.get("cp_cnvkit_scatter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_cnvkit_scatter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_cnvkit_scatter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_cnvkit_scatter", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_cnvkit_scatter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_cnvkit_scatter", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_cnvkit_scatter", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_cnvkit_scatter:
    input:
        "cnv_sv/cnvkit_scatter/{sample}_T_chr{chr}.png",
    output:
        "Results/{project}/{sample}/CNV/{sample}_T_chr{chr}.png",
    threads: config.get("cp_cnvkit_scatter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_cnvkit_scatter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_cnvkit_scatter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_cnvkit_scatter", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_cnvkit_scatter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_cnvkit_scatter", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_cnvkit_scatter", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_multiqc:
    input:
        "qc/multiqc/multiqc_DNA.html",
    output:
        "Results/MultiQC_T.html",
    threads: config.get("cp_multiqc", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_multiqc", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_multiqc", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_multiqc", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_multiqc", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_multiqc", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_multiqc", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_spring_archive:
    input:
        "compression/spring/{sample}_{flowcell}_{lane}_{barcode}_{type}.spring",
    output:
        "Archive/{project}/{sample}_{flowcell}_{lane}_{barcode}_{type}.spring",
    threads: config.get("cp_spring_archive", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_spring_archive", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_spring_archive", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_spring_archive", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_spring_archive", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_spring_archive", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_spring_archive", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"
