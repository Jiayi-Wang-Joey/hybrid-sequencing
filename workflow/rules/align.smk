configfile: "config.yaml"
align_genome = expand(res_dir + "align/genome/removed_cell/{sample}.aligned.tagged.bam", sample=SAMPLE_ALL)
############### align to transcriptome #############
rule extract_transcriptome:
    input:
        transcriptome=config["annotation"],
        genome=config["genome"]
    output:
        "data/transcriptome/transcriptome.fa"
    conda:
        "../envs/align.yaml"
    shell:
        """
        gffread -w {output} -g {input.genome} {input.transcriptome}
        """

rule add_BC_tags_to_name:
    input:
        res_dir + "preprocessing/deduplication/{sample}.unmapped.bam"
    output:
        res_dir + "align/transcriptome/keep_BC/{sample}.tags.bam"
    log:
        stderr="logs/align/transcriptome/keep_BC/{sample}.stderr"
    conda:
        "../envs/align.yaml"
    shell:
        """
        samjdk --samoutputformat BAM \
            -e '
            String cb = (String) record.getAttribute(\\"CB\\"); 
            String name = record.getReadName(); 
            record.setReadName(cb + \\"_\\" + name); 
            return record;' {input} > {output} 2> {log.stderr}
        """

rule bam2fq:
    input:
        rules.add_BC_tags_to_name.output
    output:
        res_dir + "align/transcriptome/bam2fq/{sample}.fastq.gz"
    log:
        stderr="logs/align/transcriptome/bam2fq/{sample}.stderr"
    threads: config["convert_threads"]
    params:
        threads=config["convert_threads"]
    conda:
       "../envs/align.yaml"
    shell:
        """
        samtools bam2fq -@ {params.threads} {input} | gzip > {output}
        """

rule align_minimap2_transcriptome:
    input:
        reads=rules.bam2fq.output,
        transcriptome=rules.extract_transcriptome.output,
    output:
        res_dir + "align/transcriptome/run_minimap2_novel_transcriptome/{sample}.aligned.bam",
    params:
        align_map_bam_threads=config["align_map_bam_threads"],
    threads: config["align_map_bam_threads"]
    log:
        stdout="logs/align/transcriptome/run_minimap2_transcriptome/{sample}.stdout",
        stderr="logs/align/transcriptome/run_minimap2_transcriptome/{sample}.stderr",
    conda:
        "../envs/align.yaml"
    shell:
        """
        minimap2 -ax map-hifi -N 100 --sam-hit-only --for-only \
            -t {params.align_map_bam_threads}  \
            {input.transcriptome} {input.reads} > {output}
        """

rule add_BC_tags:
    input:
        rules.align_minimap2_transcriptome.output,
    output:
        res_dir + "align/transcriptome/add_BC_tags/{sample}.aligned.tagged.bam",
    log:
        stderr="logs/align/transcriptome/add_BC_tags/{sample}.stderr",
    conda:
        "../envs/align.yaml"
    shell:
        """
        samjdk --samoutputformat BAM \
            -e '
            String s = record.getReadName(); 
            int u = s.indexOf(\\"_\\");record.setReadName(s.substring(u+1)); 
            record.setAttribute(\\"CB\\",s.substring(0,u));
            return record;' {input} > {output} 2> {log.stderr}
        """

############### align to genome #############
rule align_minimap2_genome:
    input:
        reads=rules.bam2fq.output,
        transcriptome=config["transcriptome_bed"],
        genome=config["genome"],
    output:
        res_dir + "align/genome/run_minimap2_genome/{sample}.aligned.sorted.bam"
    params:
        align_map_bam_threads=config["align_map_bam_threads"],
        align_sort_bam_threads=config["align_sort_bam_threads"],
        align_sort_bam_memory_gb=config["align_sort_bam_memory_gb"],
    threads: config["align_sort_bam_threads"] + config["align_map_bam_threads"]
    log:
        stdout="logs/align/genome/run_minimap2_genome/{sample}.stdout",
        stderr="logs/align/genome/run_minimap2_genome/{sample}.stderr"
    conda:
        "../envs/align.yaml"
    shell:
        """
        minimap2 -ax splice:hq -uf --junc-bed {input.transcriptome} \
            -t {params.align_map_bam_threads}  \
            {input.genome} {input.reads} | samtools sort \
            -@ {params.align_map_bam_threads} \
            -m{params.align_sort_bam_memory_gb}g \
            -o {output} > {log.stdout} \
            2> {log.stderr}
        """


rule add_BC_tags_genome:
    input:
        rules.align_minimap2_genome.output
    output:
        res_dir + "align/genome/add_BC_tags/{sample}.aligned.tagged.bam",
    log:
        stderr="logs/align/genome/add_BC_tags/{sample}.stderr",
    conda:
        "../envs/align.yaml"
    shell:
        """
        samjdk --samoutputformat BAM \
            -e '
            String s = record.getReadName(); 
            int u = s.indexOf(\\"_\\");record.setReadName(s.substring(u+1)); 
            record.setAttribute(\\"CB\\",s.substring(0,u));
            return record;' {input} > {output} 2> {log.stderr}
        """

def get_barcode(wc):
    base = wc.sample[:7]
    return f"{res_dir}cellranger/{base}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"


rule remove_cells:
    input:
        bam=res_dir + "align/genome/add_BC_tags/{sample}.aligned.tagged.bam",
        barcode=get_barcode
    output:
        bam = res_dir + "align/genome/removed_cell/{sample}.aligned.tagged.bam",
        bai = res_dir + "align/genome/removed_cell/{sample}.aligned.tagged.bam.bai"
    log: 
        stderr = "logs/align/genome/removed/{sample}.stderr"
    params:
        barcode = lambda wc: res_dir + f"cellranger/{wc.sample[:7]}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
    conda:
        "../envs/align.yaml"
    shell:
        """
        zcat {input.barcode} \
          | sed 's/-1$//' \
          | awk '{{print toupper($1)}}' \
          | perl -ne 'chomp; $_=reverse($_); tr/ACGTN/TGCAN/; print "CB:Z:$_\\n";' \
          > tmp/{wildcards.sample}.rc_whitelist.txt

        # Step 2: save BAM header
        samtools view -H {input.bam} > tmp/{wildcards.sample}.header.sam

        # Step 3: filter alignments by whitelist (exact CB:Z: match)
        samtools view {input.bam} \
          | LC_ALL=C grep -F -f tmp/{wildcards.sample}.rc_whitelist.txt \
          > tmp/{wildcards.sample}.body.sam

        # Step 4: combine header + body
        cat tmp/{wildcards.sample}.header.sam \
            tmp/{wildcards.sample}.body.sam \
          > tmp/{wildcards.sample}.filtered.sam

        # Step 5: convert back to BAM
        samtools view -b tmp/{wildcards.sample}.filtered.sam \
          > {output.bam} 2>> {log}

        # Step 6: index BAM
        samtools index {output.bam} 2>> {log}
        """




