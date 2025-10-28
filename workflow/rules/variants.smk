configfile: "config.yaml"
TOOL = ["Clair3-RNA", "DeepVariant"]
merged_variants = expand(res_dir + "variants/merged/{sample}.pileup.vcf.gz", sample=SAMPLE_ALL)
clinvar = expand(res_dir + "variants/clinvar/{sample}.clinvar.txt", sample=SAMPLE_ALL)
# Germline variant calling ==========================================
rule run_clair3_rna:
    input:
        bam=res_dir + "align/genome/removed_cell/{sample}.aligned.tagged.bam",
        ref=config["genome"],
        img=config["clair3-rna"]
    output:
        res_dir + "variants/Clair3-RNA/germline/filter0/{sample}.vcf.gz"
    log: 
        "logs/run_clair3_rna_{sample}.log"
    params:
        threads = 10,
        dir = directory(res_dir + "variants/Clair3-RNA/germline/{sample}")
    shell:
        """
        samtools index {input.bam}
        singularity exec -B {input.bam},{input.ref} {input.img} \
            /bin/bash -c "source /opt/conda/bin/activate /opt/conda/envs/clair3_rna && \
            /opt/bin/run_clair3_rna \
            --bam_fn {input.bam} \
            --ref_fn {input.ref} \
            --threads {params.threads} \
            --platform hifi_mas_minimap2 \
            --tag_variant_using_readiportal \
            --enable_phasing_model \
            --output_dir {params.dir} \
            --conda_prefix /opt/conda/envs/clair3_rna" \
        > {log} 2>&1
        mv {params.dir}/output.vcf.gz {output}
        """

rule run_deep_variant:
    input: 
        bam = res_dir + "align/genome/removed_cell/{sample}.aligned.tagged.bam",
        ref = config["genome"],
        img = config["deepvariant"]
    output:
        vcf = res_dir + "variants/DeepVariant/germline/filter0/{sample}.vcf.gz"
    log:
        "logs/run_deep_variant_{sample}.log"
    params: 
        threads = 8,
        tmp_dir = "/home/jiayiwang/tmp/{sample}",
        pseudoq_bam = lambda wildcards, input: input.bam.replace(".tagged.bam", ".pseudoQ.bam")
    conda:
        "../envs/variants.yaml"
    shell:
        r"""
        # 1. Create pseudo-quality BAM
        python3 - <<'EOF'
import pysam
in_bam  = "{input.bam}"
out_bam = "{params.pseudoq_bam}"
with pysam.AlignmentFile(in_bam, "rb") as infile, \
     pysam.AlignmentFile(out_bam, "wb", template=infile) as outfile:
    for read in infile:
        if read.query_length > 0:
            read.query_qualities = [93] * read.query_length
        outfile.write(read)
print("Pseudo-quality BAM created:", out_bam)
EOF

        # 2. Index the pseudoQ BAM
        samtools index {params.pseudoq_bam}

        # 3. Run DeepVariant
        singularity exec --bind /usr/lib/locale/ {input.img} /opt/deepvariant/bin/run_deepvariant \
            --model_type MASSEQ \
            --ref {input.ref} \
            --reads {params.pseudoq_bam} \
            --output_vcf {output.vcf} \
            --num_shards {params.threads} \
            --intermediate_results_dir {params.tmp_dir} \
        > {log} 2>&1
        """


# Variants filtering ==========================================
rule fil1_cov:
    input: 
        res_dir + "variants/{tool}/germline/filter0/{sample}.vcf.gz"
    output: 
        res_dir + "variants/{tool}/germline/filter1/{sample}.vcf.gz"
    params:
        dp=5
    conda:
        "../envs/variants.yaml"
    log:    
        "logs/fil_coverage_{tool},{sample}.log"
    shell:
        '''
        bcftools view -i 'DP>={params.dp}' {input} -Oz -o {output}
        '''


rule fil2_dat:
    input:
        vcf=rules.fil1_cov.output,  
        dbSNP="/home/jiayiwang/variants/data/dbSNP/00-common_all.vcf.gz",
    output:
        vcf=res_dir + "variants/{tool}/germline/filter2/{sample}.vcf.gz"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/filter_dat_{tool}_{sample}.log"
    shell:
        """
        # Index both files to ensure isec works
        bcftools index -f {input.vcf}

        # Filter out variants found in dbSNP
        bcftools isec -C {input.vcf} {input.dbSNP} -Oz -p isec_{wildcards.sample}

        # Compress and output the filtered VCF
        bcftools view -Oz -o {output.vcf} isec_{wildcards.sample}/0000.vcf.gz
        
        # Clean up
        rm -rf isec_{wildcards.sample}
        """


rule fil3_pass:
    input:
        rules.fil2_dat.output
    output:
        res_dir + "variants/{tool}/germline/filter3/{sample}.vcf.gz"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/filter_pass_{tool},{sample}.log"
    shell:
        """
        bcftools index -f {input}
        bcftools view -i 'FILTER="PASS" || QUAL>=10' {input} -Oz -o {output}
        bcftools index -f {output}
        """


rule fil4_cds:
    input:
        vcf=rules.fil3_pass.output,
        bed=config["cds_bed"]
    output:
        res_dir + "variants/{tool}/germline/filter4/{sample}.vcf"
    log:
        "logs/filter_cds_{tool},{sample}.log"
    conda:
        "../envs/variants.yaml"
    shell:
        """
        bcftools index -f {input.vcf}
        bcftools view -R {input.bed} {input.vcf} -o {output}
        """
    
rule anno_vcf:
    input:
        vcf = rules.fil4_cds.output
    output:
        res_dir + "variants/{tool}/germline/filter4/{sample}.anno.vcf.gz"
    log:
        "logs/anno_vcf_{tool},{sample}.log"
    params:
        expr = " | ".join([
            f"(ANN[*].GENE = '{g}')" for g in [
                "ATP1A1","ATP2B3","CACNA1D","CACNA1H","KCNJ5","GNA11","MCOLN3","CLCN2","CTNNB1","CYP11B2","FRRS1L",
                "CELA1","UBE2Q2L","UBTD1","NCAM1","CAMTA2","TMEM67","KCNA4","KCNH3","KCNK5","PER1","WNK2","KCNMB3","INSR",
                "PDE8A","PPP1R9A","USP3","NR0B2","AQP11","NFATC4","STAR","CYP11B1","KCNK3","SLC8A1","CALM1","CREB1","PRKACA","MTOR",
                "GNAQ","SLC24A1","PEX1","VDR","TSC2","CDKN2B","MAP3K6","SLC30A1","NR5A1","HSD3B2","NR4A2","WNK1","MCLON3",
            ]
        ])
    shell:
        """
        java -Xmx8g -jar ~/snpEff/snpEff.jar GRCh38.86 {input.vcf} \
        | java -jar ~/snpEff/SnpSift.jar filter "{params.expr}" \
        | bgzip -c > {output}
        """



rule merge_vcf:
    input:
        lambda wildcards: expand(
            res_dir + "variants/{tool}/germline/filter4/{sample}.anno.vcf.gz",
            tool=TOOL, sample=wildcards.sample
        )
    output:
        res_dir + "variants/merged/{sample}.anno.vcf.gz"
    log:
        "logs/variants/merged/{sample}.merge.log"
    conda:
        "../envs/variants.yaml"
    threads: 4
    shell:
        """
        for f in {input}; do
            bcftools index --force "$f" >> {log} 2>&1
        done
        bcftools merge -Oz -o {output} {input} >> {log} 2>&1
        bcftools index --force {output} >> {log} 2>&1
        """

rule vcf_pileup_annotate:
    input:
        vcf = res_dir + "variants/merged/{sample}.anno.vcf.gz",
        bam = res_dir + "align/genome/removed_cell/{sample}.aligned.tagged.bam",
        ref = config["genome"]  
    output:
        vcf = res_dir + "variants/merged/{sample}.pileup.vcf.gz"
    log:
        "logs/pileup/{sample}.log"
    threads: 8
    conda:
        "../envs/variants.yaml"
    shell:
        """
        # Index BAM and input VCF if missing
        [ -f {input.bam}.bai ] || samtools index {input.bam}
        [ -f {input.vcf}.tbi ] || bcftools index -f {input.vcf}

        mkdir -p tmp

        # Extract variant positions from input VCF
        bcftools query -f '%CHROM\t%POS0\t%POS\n' {input.vcf} > tmp/{wildcards.sample}.sites.bed

        # Step 1: mpileup BAM at variant positions → plain VCF (adds DP/AD)
        bcftools mpileup \
            -f {input.ref} \
            -a FORMAT/DP,FORMAT/AD \
            -R tmp/{wildcards.sample}.sites.bed \
            -Ov {input.bam} \
            -o tmp/{wildcards.sample}.pileup.vcf 2>> {log}

        # Step 2: bgzip and index temporary VCF
        bgzip -c tmp/{wildcards.sample}.pileup.vcf > tmp/{wildcards.sample}.pileup.vcf.gz
        tabix -f -p vcf tmp/{wildcards.sample}.pileup.vcf.gz

        # Step 3: annotate original ALT alleles from input VCF
        bcftools annotate \
            -a {input.vcf} \
            -c CHROM,POS,REF,ALT \
            tmp/{wildcards.sample}.pileup.vcf.gz \
            -Oz -o tmp/{wildcards.sample}.pileup.with_ALT.vcf.gz 2>> {log}

        # Index intermediate annotated VCF
        bcftools index -f tmp/{wildcards.sample}.pileup.with_ALT.vcf.gz

        # Step 4: SnpEff annotation (canonical transcripts only)
        java -Xmx8g -jar ~/snpEff/snpEff.jar GRCh38.86 -canon tmp/{wildcards.sample}.pileup.with_ALT.vcf.gz \
        | bgzip -c > {output.vcf} 2>> {log}

        # Index final annotated VCF
        bcftools index -f {output.vcf}

        # Cleanup temporary files
        rm tmp/{wildcards.sample}.sites.bed \
           tmp/{wildcards.sample}.pileup.vcf \
           tmp/{wildcards.sample}.pileup.vcf.gz \
           tmp/{wildcards.sample}.pileup.vcf.gz.tbi \
           tmp/{wildcards.sample}.pileup.with_ALT.vcf.gz \
           tmp/{wildcards.sample}.pileup.with_ALT.vcf.gz.csi
        """



rule clinvar:
    input:
        vcf=res_dir + "variants/merged/{sample}.pileup.vcf.gz",
        ann="data/annotation/clinvar.vcf.gz"
    output:
        vcf=res_dir + "variants/clinvar/{sample}.clinvar.vcf",
        txt=res_dir + "variants/clinvar/{sample}.clinvar.txt"
    conda:
        "../envs/variants.yaml"
    log:
        "logs/clinvar_{sample}.log"
    shell:
        """
        java -jar ~/snpEff/SnpSift.jar annotate {input.ann} {input.vcf} > {output.vcf}
        echo -e "CHROM\tPOS\tREF\tALT\tTYPE\tID\tALLELEID\tCLNDN\tCLNSIG\tCLNSIGCONF\tCLNSIGINCL\tCLNVC\tGENEINFO\tEFFECT\tAD\tGQ\tGT\tVAF\tDP\tVCF_ID" > {output.txt}

        # Conditional check for tool type

        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%ID\t%INFO/ALLELEID\t%INFO/CLNDN\t%INFO/CLNSIG\t%INFO/CLNSIGCONF\t%INFO/CLNSIGINCL\t%INFO/CLNVC\t%INFO/GENEINFO\t%INFO/ANN[*].EFFECT\t[%AD\t%AD\t%DP]' {output.vcf} >> {output.txt}
    

        awk 'BEGIN {{FS="\\t"; OFS="\\t"}} {{if (NR == 1) print $0; else if ($6 != ".") print $0, NR - 1}}' {output.txt} > {output.txt}.tmp && mv {output.txt}.tmp {output.txt}
        """







