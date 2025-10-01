configfile: "config.yaml"

dedup = expand(res_dir + "preprocessing/deduplication/{sample}_all.unmapped.bam", sample=SAMPLE_LR) 


rule primer_removal:
    priority: 99
    input:
        reads=config["lr_bam_dir"] +  "{sample}-segmented.bam",
        primers=config["primers"]
    output:
        res_dir + "preprocessing/primer_removal/{sample}.fl.5p--3p.bam"
    log:
        stderr="logs/preprocessing/primer_removal/{sample}.stderr"
    params:
        threads=config["lima_threads"],
        filename=res_dir + "preprocessing/primer_removal/{sample}.fl.bam"
    conda:
        "../envs/preprocessing.yaml"
    shell:
        """
        lima {input.reads} {input.primers} {params.filename} --no-reports --num-threads {params.threads} --isoseq
        """

rule tag:
    priority: 98
    input:
        rules.primer_removal.output
    output:
        res_dir + "preprocessing/tags/{sample}.flt.bam"
    log:
        stderr="logs/preprocessing/tags/{sample}.stderr"
    params:
        threads=config["tag_threads"],
    conda:
        "../envs/preprocessing.yaml"
    shell:
        """
        isoseq tag {input} {output} --design T-12U-16B --num-threads {params.threads}
        """

rule refine:
    priority: 97
    input:
        reads=rules.tag.output,
        primers=config["primers"],
    output:
        res_dir + "preprocessing/refine/{sample}.fltnc.bam",
    log:
        stderr="logs/preprocessing/refine/{sample}.stderr",
    params: 
        threads=config["refine_threads"]
    conda:
        "../envs/preprocessing.yaml"
    shell:
        """
        isoseq refine {input.reads} {input.primers} {output} --require-polya --num-threads {params.threads}
        """


rule BC_correction:
    priority: 96
    input:
        reads=rules.refine.output
    output:
        res_dir + "preprocessing/BC_correction/{sample}.corrected.bam",
    log:
        stderr="logs/preprocessing/BC_correction/{sample}.stderr"
        
    params: 
        threads = config["bc_corr_threads"],
        barcodes = config["whitelist"]

    conda:
        "../envs/preprocessing.yaml"
    shell:
        """
        isoseq correct  --method percentile --percentile 95 --barcodes {params.barcodes} \
            --num-threads {params.threads} {input.reads} {output}
        """

rule sort_bam:
    priority: 95
    input:
        rules.BC_correction.output,
    output:
        res_dir + "preprocessing/sorting/{sample}.sorted.bam",
    log:
        stderr="logs/preprocessing/sorting/{sample}.stderr",
    params:
        threads=config["sorting_threads"],
        memory=config["sorting_memory"],
    conda:
        "../envs/preprocessing.yaml"
    shell:
        """
        samtools sort -t CB -@ {params.threads} -m{params.memory}g {input} -o {output}
        """

rule deduplication:
    priority: 94
    input:
        res_dir + "preprocessing/sorting/{sample}.sorted.bam",
    output:
        res_dir + "preprocessing/deduplication/{sample}.unmapped.bam"
    log:
        stderr="logs/preprocessing/deduplication/{sample}.stderr"
    params:
        threads=config["dedup_threads"],
    conda:
        "../envs/preprocessing.yaml"
    shell:
        """
        isoseq groupdedup --num-threads {params.threads} {input} {output}
        """

rule deduplication_all:
    priority: 94
    input:
        res_dir + "preprocessing/sorting/{sample}.sorted.bam",
    output:
        res_dir + "preprocessing/deduplication/{sample}_all.unmapped.bam"
    log:
        stderr="logs/preprocessing/deduplication/{sample}.stderr"
    params:
        threads=config["dedup_threads"],
    conda:
        "../envs/preprocessing.yaml"
    shell:
        """
        isoseq groupdedup --keep-non-real-cells --num-threads {params.threads} {input} {output}
        """