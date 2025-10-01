configfile: "config.yaml"

# qualimap = expand(
#     "results/read_qc/qualimap_{sample}_lr/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt",
#     sample=SAMPLE
# ) + expand(
#     "results/read_qc/qualimap_{sample}_sr/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt",
#     sample=SAMPLE
# )

read_distribution = expand("results/read_qc/read_distribution_{sample}_lr.tsv", sample=SAMPLE) + \
    expand("results/read_qc/read_distribution_{sample}_sr.tsv", sample=SAMPLE)

rule plot_gene_coverage:
    input:
        "workflow/scripts/plt-gene_coverage.R",
    output:
        "plts/F1/gene_coverage.pdf"
    params:
        lambda wc, input: ";".join(genebody_coverage)
    log:
        "logs/plt-gene_coverage.log"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        {params} {output[0]}" {input[0]} {log}'''

rule plot_read_distribution:
    input:
        "workflow/scripts/plt-read_distribution.R",
    output:
        "plts/F1/read_distribution.pdf"
    params:
        lambda wc, input: ";".join(read_distribution)
    log:
        "logs/plt-read_distribution.log"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        {params} {output[0]}" {input[0]} {log}'''