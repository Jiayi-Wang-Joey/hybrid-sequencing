
SEED = [42, 2026, 2, 7025, 226]
SAMPLE = ["sampleA", "sampleB"]
downsample_fastq = expand("/home/jiayiwang/loggedfs_david_penton/result/downsample/fastq/{sample}_{seed}_R1.fastq.gz", sample=SAMPLE, seed=SEED) +\
    expand("/home/jiayiwang/loggedfs_david_penton/result/downsample/fastq/{sample}_{seed}_R2.fastq.gz", sample=SAMPLE, seed=SEED)
downsample_count =expand(
            "results/cellranger/{sample}_{seed}",
            seed=SEED,
            sample=SAMPLE
        )

rule downsample_fastq:
    input:
        R1=lambda wc: {
            "sampleA": "/home/jiayiwang/loggedfs_david_penton/data/NovaSeq/387241_2-patient4_normal/387241_2-patient4_normal_S4_L008_R1_001.fastq.gz",
            "sampleB": "/home/jiayiwang/loggedfs_david_penton/data/NovaSeq/387241_1-patient4_tumor/387241_1-patient4_tumor_S1_L008_R1_001.fastq.gz"
        }[wc.sample],
        R2=lambda wc: {
            "sampleA": "/home/jiayiwang/loggedfs_david_penton/data/NovaSeq/387241_2-patient4_normal/387241_2-patient4_normal_S4_L008_R2_001.fastq.gz",
            "sampleB": "/home/jiayiwang/loggedfs_david_penton/data/NovaSeq/387241_1-patient4_tumor/387241_1-patient4_tumor_S1_L008_R2_001.fastq.gz"
        }[wc.sample]
    output:
        R1="/home/jiayiwang/loggedfs_david_penton/result/downsample/fastq/{sample}_{seed}/{sample}_{seed}_S1_L001_R1_001.fastq.gz",
        R2="/home/jiayiwang/loggedfs_david_penton/result/downsample/fastq/{sample}_{seed}/{sample}_{seed}_S1_L001_R2_001.fastq.gz"
    params:
        n=lambda wc: {
            "sampleA": 85182821,
            "sampleB": 78196246
        }[wc.sample]
    log:
        "logs/downsample_fastq/{sample}_{seed}.log"
    shell:
        """
        seqtk sample -s {wildcards.seed} {input.R1} {params.n} | gzip > {output.R1} 2> {log}
        seqtk sample -s {wildcards.seed} {input.R2} {params.n} | gzip > {output.R2} 2> {log}
        """


rule cellranger_count_downsample:
    input:
        R1 = "/home/jiayiwang/loggedfs_david_penton/result/downsample/fastq/{sample}_{seed}/{sample}_{seed}_S1_L001_R1_001.fastq.gz",
        R2 = "/home/jiayiwang/loggedfs_david_penton/result/downsample/fastq/{sample}_{seed}/{sample}_{seed}_S1_L001_R2_001.fastq.gz"
    output:
        directory("results/cellranger/{sample}_{seed}")
    params:
        cores = 8,
        fastq_dir = "/home/jiayiwang/loggedfs_david_penton/result/downsample/fastq/{sample}_{seed}",
        cellranger = config["cellranger"],
        transcriptome = config["cr_transcriptome"]
    log:
        "logs/cellranger/{sample}_{seed}.log"
    shell:
        """
        {params.cellranger} count \
          --id={wildcards.sample}_{wildcards.seed} \
          --transcriptome={params.transcriptome} \
          --output-dir={output} \
          --fastqs={params.fastq_dir} \
          --sample={wildcards.sample}_{wildcards.seed} \
          --localcores={params.cores} \
          --localmem=64 \
          --create-bam=false \
          2> {log}
        """