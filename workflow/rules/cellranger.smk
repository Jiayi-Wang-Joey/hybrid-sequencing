rule run_cellranger:
    input:
        config["sr_fastq_dir"] + "{sample}"
    output:
        directory(res_dir + "cellranger/{sample}")
    conda:
        "envs/cellranger.yaml"
    params:
        cores = 8
    log:
        "logs/cellranger/{sample}.log"
    shell:
        """
        cellranger count \
          --id={wildcards.sample} \
          --transcriptome=GRCh38p14_gencode_v47 \
          --fastqs={input} \
          --sample={wildcards.sample} \
          --localcores={params.cores} \
          --localmem=64 \
          --create-bam=true \
          2> {log}
        """


rule removed_cell:
    input:
        bam = res_dir + "cellranger/{sample}/outs/possorted_genome_bam.bam",
        barcode = res_dir + "cellranger/{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
    output:
        bam = res_dir + "cellranger/{sample}/outs/filtered_genome_bam.bam"
    conda:
        "../envs/samtools.yaml"
    log:
        "logs/cellranger/removecells.{sample}.log"
    shell:
        """
        # Step 1: extract all barcodes present in BAM (CB tag)
        samtools view {input.bam} \
          | awk '{{ for (i=12;i<=NF;i++) if ($i ~ /^CB:Z:/) print substr($i,6) }}' \
          | sort -u > tmp/{wildcards.sample}.all_barcodes.txt

        # Step 2: prepare whitelist of barcodes to keep
        zcat {input.barcode} > tmp/{wildcards.sample}.rc_whitelist.txt

        # Step 3: save BAM header
        samtools view -H {input.bam} > tmp/{wildcards.sample}.header.sam

        # Step 4: filter alignments by whitelist (exact CB:Z: match)
        samtools view {input.bam} \
          | LC_ALL=C grep -F -f tmp/{wildcards.sample}.rc_whitelist.txt \
          > tmp/{wildcards.sample}.body.sam

        # Step 5: combine header + body
        cat tmp/{wildcards.sample}.header.sam \
            tmp/{wildcards.sample}.body.sam \
          > tmp/{wildcards.sample}.filtered.sam

        # Step 6: convert back to BAM
        samtools view -b tmp/{wildcards.sample}.filtered.sam \
          > {output.bam} 2>> {log}

        # Step 7: index BAM
        samtools index {output.bam} 2>> {log}

        # Step 8: cleanup
        rm -f tmp/{wildcards.sample}.all_barcodes.txt \
              tmp/{wildcards.sample}.rc_whitelist.txt \
              tmp/{wildcards.sample}.header.sam \
              tmp/{wildcards.sample}.body.sam \
              tmp/{wildcards.sample}.filtered.sam
        """
