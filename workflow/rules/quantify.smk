configfile: "config.yaml"
count = expand(res_dir+"quantify/transcriptome/oarfish/{sample}.count.mtx", sample=SAMPLE_ALL)
sce = expand(res_dir+"quantify/transcriptome/gene_sce/{sample}.rds", sample=SAMPLE_ALL)
rule sort_bams_by_barcode_transcriptome:
    input:
        res_dir + "align/transcriptome/add_BC_tags/{sample}.aligned.tagged.bam",
    output:
        res_dir + "quantify/transcriptome/sort_bams_by_barcode/{sample}.aligned.sorted_barcode.bam",
    threads: config["sorting_threads"]
    params:
        sorting_threads=config["sorting_threads"],
        sorting_memory=config["sorting_memory"],
    log:
        "logs/quantify/transcriptome/sort_bams_by_barcode/{sample}.out",
    conda:
        "../envs/quantify.yaml"
    shell:
        """
        samtools sort -@ {params.sorting_threads} -m{params.sorting_memory}g -t CB {input} -o {output} &> {log}
        """

rule oarfish_quantify:
    input:
        res_dir+"quantify/transcriptome/sort_bams_by_barcode/{sample}.aligned.sorted_barcode.bam"
    output:
        res_dir+"quantify/transcriptome/oarfish/{sample}.count.mtx",
        res_dir+"quantify/transcriptome/oarfish/{sample}.barcodes.txt",
        res_dir+"quantify/transcriptome/oarfish/{sample}.features.txt"
    log:
        "logs/quantify/transcriptome/oarfish_quantify/{sample}.out"
    conda:
        "../envs/quantify.yaml"
    params:
        outdir = res_dir + "quantify/transcriptome/oarfish/{sample}",
        thread=20
    shell:
        """
        oarfish --single-cell --model-coverage --filter-group no-filters --output {params.outdir} --alignments {input} --threads {params.thread} &> {log}
        """

rule quantify_extract_transcript_gene_mapping:
    input:
        config["annotation"]
    output:
        "data/extract_transcript_gene_mapping/transcript_gene_mapping.tsv",
    threads: 1
    log:
        "logs/quantify/transcriptome/extract_transcript_gene_mapping/log.out",
    conda:
        "../envs/quantify.yaml"
    shell:
        """
        gffread {input} --table transcript_id,gene_id > {output} 2> {log};
        sed -i '1i \ttranscript_id\tgene_id' {output} 2>> {log}
        """

rule generate_sce:
    input:
        "workflow/scripts/generate_sce.R",
        tx="data/extract_transcript_gene_mapping/transcript_gene_mapping.tsv",
        mtx=res_dir+"quantify/transcriptome/oarfish/{sample}.count.mtx",
        fts=res_dir+"quantify/transcriptome/oarfish/{sample}.features.txt",
        bcd=res_dir+"quantify/transcriptome/oarfish/{sample}.barcodes.txt"
    output:
        tsce=res_dir+"quantify/transcriptome/transcript_sce/{sample}.rds",
        gsce=res_dir+"quantify/transcriptome/gene_sce/{sample}.rds"
    log:
        "logs/quantify/transcriptome/generate_sce/{sample}.out"
    shell:
        """
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        mtx={input.mtx} fts={input.fts} bcd={input.bcd} tx={input.tx} tsce={output.tsce} gsce={output.gsce}" {input[0]} {log}
        """