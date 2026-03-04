import json
import itertools

configfile: "config.yaml"
R = config["R"]
res_dir = config["res_dir"]

# wildcards
SAMPLE = ["sampleA", "sampleB"]
SAMPLE_LR = ["sampleA", "sampleB", "sampleA_twist", "sampleB_twist"]
SAMPLE_ALL = ["sampleA_twist_all", "sampleB_twist_all", "sampleA_all", "sampleB_all"]

#rules
include: "workflow/rules/cellranger.smk"
include: "workflow/rules/preprocess.smk"
include: "workflow/rules/align.smk"
include: "workflow/rules/quantify.smk"
include: "workflow/rules/variants.smk"
#include: "workflow/rules/downsample.smk"


# results 
################# preprocess ####################
preprocess = {
    "dedup": dedup
}
################## Alignment #####################
align = {
    "align_genome": align_genome
}

################## Quantification ################
quantify = {
    "sce": sce
    #"sce_raw": sce_raw
}

################## Variant calling ##################
all_variants = expand(res_dir + "variants/merged/{sample}.pileup.vcf.gz", sample=SAMPLE_LR)
variants = {
    "merged_variants": merged_variants,
    "clinvar": clinvar
}

################ Revision ##############
# downsample = {
#     "downsample_fastq": downsample_fastq,
#     "downsample_count": downsample_count
# }

############## set up ##################
rule all:
    input:
        [x for x in preprocess.values()],
        [x for x in align.values()],
        [x for x in quantify.values()],
        [x for x in variants.values()]
        #[x for x in downsample.values()],



