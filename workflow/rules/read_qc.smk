configfile: "config.yaml"
genebody_coverage = expand("results/read_qc/{sample}_sr.geneBodyCoverage.txt", sample=SAMPLE) + \
    expand("results/read_qc/{sample}_lr.geneBodyCoverage.txt", sample=SAMPLE)


rule subsample_gtf:
    input:
        "workflow/scripts/subsample_transcripts.R",
        config["annotation"]
    output:
        "data/transcriptome/subsample_transcripts.gtf"
    params:
        n=10000
    log:
        "logs/read_qc/subsample_gtf.log"
    shell:
        '''
        {R} CMD BATCH --no-restore --no-save "--args\
        {params} res={output[0]} gtf={input[1]}" {input[0]} {log}'''

rule convert_gtf_to_bed:
    input:
        rules.subsample_gtf.output
    output:
        "data/transcriptome/transcriptome.bed"
    log:
        "logs/read_qc/convert_gtf_to_bed.log"
    conda:
        "../envs/read_qc.yaml"
    shell:
        """
        cat {input} 2>> {log} |\
            gtfToGenePred /dev/stdin /dev/stdout 2>> {log} |\
            genePredToBed /dev/stdin /dev/stdout > \
            {output} 2>> {log}
        """
        
# rule subset_bed:
#     input:
#         rules.convert_gtf_to_bed.output
#     output:
#         "results/read_qc/gencode.v47.subset10000.bed"
#     log:
#         "logs/subset_bed.log"
#     shell:
#         """
#         shuf -n 10000 {input} > {output} 2> {log}
#         """



rule read_distribution_lr:
    input:
        bam = res_dir + "align/genome/add_BC_tags/{sample}.aligned.tagged.bam",
        bed = "/home/jiayiwang/tools/RSeQC-5.0.1/gene_model/hg38_GENCODE_V47.bed"
    output:
        "results/read_qc/read_distribution_{sample}_lr.tsv"
    conda:
        "../envs/read_qc.yaml"
    log:
        "logs/read_distribution_{sample}_lr.log"
    shell:
        """
        read_distribution.py -i {input.bam} -r {input.bed} > {output} 2> {log}
        """

rule read_distribution_sr:
    input:
        bam = res_dir + "cellranger/{sample}/outs/possorted_genome_bam.bam",
        bed = "/home/jiayiwang/tools/RSeQC-5.0.1/gene_model/hg38_GENCODE_V47.bed"
    output:
        "results/read_qc/read_distribution_{sample}_sr.tsv"
    conda:
        "../envs/read_qc.yaml"
    log:
        "logs/read_distribution_{sample}_sr.log"
    shell:
        """
        read_distribution.py -i {input.bam} -r {input.bed} > {output} 2> {log}
        """

rule gene_coverage_lr:
    input:
         bam = res_dir + "align/genome/removed_cell/{sample}_all.aligned.tagged.bam",
         bed = rules.convert_gtf_to_bed.output
    output:
        "results/read_qc/{sample}_lr.geneBodyCoverage.txt"
    params:
        "results/read_qc/{sample}_lr"
    conda:
        "../envs/read_qc.yaml"
    log:
        "logs/geneBody_coverage/lr_{sample}.log"
    shell:
        """
        geneBody_coverage.py -i {input.bam} -r {input.bed} -o {params} 2> {log}
        """


rule gene_coverage_sr:
    input:
         bam = res_dir + "cellranger/{sample}/outs/filtered_genome_bam.bam",
         bed = rules.convert_gtf_to_bed.output
    output:
        "results/read_qc/{sample}_sr.geneBodyCoverage.txt"
    params:
        "results/read_qc/{sample}_sr"
    conda:
        "../envs/read_qc.yaml"
    log:
        "logs/geneBody_coverage/sr_{sample}.log"
    shell:
        """
        geneBody_coverage.py -i {input.bam} -r {input.bed} -o {params} 2> {log}
        """

# rule qualimap_lr:
#     input:
#         bam = res_dir + "align/genome/add_BC_tags/{sample}.aligned.tagged.bam",
#         gtf = config["annotation"]
#     output:
#         "results/read_qc/qualimap_{sample}_lr/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt"
#     log:
#         "logs/qualimap_{sample}_sr.log"
#     params:
#         threads = 8,
#         dir = "results/read_qc/qualimap_{sample}_lr/"
#     conda:
#         "../envs/read_qc.yaml"
#     shell:
#         """
#         if [ ! -f {input.bam}.bai ]; then
#             samtools index {input.bam} >> {log} 2>&1
#         fi
#         qualimap rnaseq  \
#         -a proportional \
#         -bam {input.bam} \
#         -gtf {input.gtf} \
#         -outdir {params.dir} \
#         --java-mem-size=16G \
#         -oc {params.threads} \
#         2>> {log}
#         """

# rule qualimap_sr:
#     input:
#         bam = res_dir + "cellranger/{sample}/outs/possorted_genome_bam.bam",
#         gtf = config["annotation"]
#     output:
#         "results/read_qc/qualimap_{sample}_sr/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt"
#     log:
#         "logs/qualimap_{sample}_sr.log"
#     params:
#         threads = 8,
#         dir = "results/read_qc/qualimap_{sample}_sr/"
#     conda:
#         "../envs/read_qc.yaml"
#     shell:
#         """
#         if [ ! -f {input.bam}.bai ]; then
#             samtools index {input.bam} >> {log} 2>&1
#         fi
#         qualimap rnaseq  \
#         -a proportional \
#         -bam {input.bam} \
#         -gtf {input.gtf} \
#         -outdir {params.dir} \
#         --java-mem-size=16G \
#         -oc {params.threads} \
#         2>> {log}
#         """


rule barcode_count_lr:
    input:
        res_dir + "preprocessing/sorting/{sample}.sorted.bam"
    output:
        "results/read_qc/barcode_count_{sample}_lr.txt"
    log:
        "logs/read_qc/barcode_count_{sample}.log"
    conda:
        "../envs/read_qc.yaml"
    shell:
        """
        samtools view {input} 2> {log} \
        | awk '{{for(i=12;i<=NF;i++){{if($i ~ /^CB:Z:/){{print substr($i,6)}}}}}}' 2>> {log} \
        | sort 2>> {log} \
        | uniq -c 2>> {log} \
        > {output}
        """

rule barcode_count_sr:
    input:
        fastq_r1 = config["sr_fastq_dir"] + "{sample}_R1.fastq.gz"
    output:
        "results/read_qc/barcode_count_{sample}_sr.txt"
    log:
        "logs/read_qc/barcode_count_{sample}_lr.log"
    conda:
        "../envs/read_qc.yaml"
    shell:
        """
        zcat {input.fastq_r1} 2>> {log} \
        | awk 'NR % 4 == 2 {{print substr($0,1,16)}}' 2>> {log} \
        | sort 2>> {log} \
        | uniq -c 2>> {log} \
        > {output}
        """

