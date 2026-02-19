import json
import itertools
configfile: "config.yaml"
R = config["R"]

# WILDCARDS --------------------------------------------------------------------
ISOSEQ_SAMPLE = ["HG004-Baylor-IsoSeq",
   "HG005-Baylor-IsoSeq", "HG002-Baylor-IsoSeq"
]
MASSEQ_SAMPLE = ["HG004-MasSeq",
"HG005-MasSeq","HG002-MasSeq"
]
PACBIO_SAMPLE = ISOSEQ_SAMPLE + MASSEQ_SAMPLE
ONT_SAMPLE = [
    "HG004-cDNAxR09", "HG004-cDNAxR10",  "HG004-dRNA004" ,"HG004-dRNA002",
    "HG005-cDNAxR09", "HG005-cDNAxR10",  "HG005-dRNA004" ,"HG005-dRNA002",
    "HG002-cDNAxR10", "HG002-dRNA004", "HG002-dRNA002", "HG002-cDNAxR09"
    ] 

BAM = ["origin"] # transformation of bam
PACBIO_ALIGNER = ["minimap2", "pbmm2"] # pacbio aligner
ONT_ALIGNER = ["minimap2"] # ONT aligner
COVERAGE = [1, 5, 10, 30, 50, 100] # coverage
TOOL = ["Clair3-RNA", "DeepVariant", "longcallR", "longcallR-nn", "GATK"] # variant caller
PHASER = ["WhatsHap", "longphase", "HapCUT2"]

ISOSEQ_STATUS = ["raw", "preprocessed"]

# RULES --------------------------------------------------------------------
include: "workflow/rules/align.smk"
include: "workflow/rules/variant_calling.smk"
include: "workflow/rules/variant_benchmark.smk"
include: "workflow/rules/phasing.smk"


# RESULTS --------------------------------------------------------------------

#### Alignment 
preprocess = expand("data/preprocessed/{sample}.bam", sample=ISOSEQ_SAMPLE)

origin = (expand("results/align/{status}/origin/{aligner}_{sample}.aligned.bam",
                aligner=PACBIO_ALIGNER, sample=ISOSEQ_SAMPLE, status=ISOSEQ_STATUS) + 
         expand("results/align/raw/origin/{aligner}_{sample}.aligned.bam",
                aligner=PACBIO_ALIGNER, sample=MASSEQ_SAMPLE) + 
         expand("results/align/raw/origin/{aligner}_{sample}.aligned.bam",
                aligner=ONT_ALIGNER, sample=ONT_SAMPLE))

transformed = expand("results/align/{status}/transformed/{aligner}_{sample}.aligned.bam", 
    aligner=PACBIO_ALIGNER, sample=ISOSEQ_SAMPLE, status=ISOSEQ_STATUS) + \
    expand("results/align/raw/transformed/{aligner}_{sample}.aligned.bam",
                aligner=PACBIO_ALIGNER, sample=MASSEQ_SAMPLE) + \
    expand("results/align/raw/transformed/{aligner}_{sample}.aligned.bam", aligner=ONT_ALIGNER, sample=ONT_SAMPLE)

qual = expand("results/align/{status}/quality/{aligner}_{sample}.tsv", aligner=PACBIO_ALIGNER, sample=ISOSEQ_SAMPLE, status=ISOSEQ_STATUS) + \
    expand("results/align/raw/quality/{aligner}_{sample}.tsv", aligner=PACBIO_ALIGNER, sample=MASSEQ_SAMPLE) + \
    expand("results/align/raw/quality/{aligner}_{sample}.tsv", aligner=ONT_ALIGNER, sample=ONT_SAMPLE) 

happy = expand("results/happy/{tool}_origin_{aligner}_{sample}_{coverage}.summary.csv", tool=TOOL, aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, coverage=COVERAGE) + \
   expand("results/happy/{tool}_origin_{aligner}_{sample}_{coverage}.summary.csv", tool=TOOL,  aligner=ONT_ALIGNER, sample=ONT_SAMPLE, coverage=COVERAGE)

happy_transformed = expand("results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.summary.csv", tool=TOOL, aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, coverage=COVERAGE, bamtype=BAM) + \
   expand("results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.summary.csv", tool=TOOL,  aligner=ONT_ALIGNER, sample=ONT_SAMPLE, coverage=COVERAGE, bamtype=BAM)

# roc = expand("results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.roc.csv.gz",  tool=TOOL, aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, bamtype=BAM,  coverage=COVERAGE) + \
#     expand("results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.roc.csv.gz",  tool=TOOL, aligner=ONT_ALIGNER, sample=ONT_SAMPLE, bamtype=BAM,  coverage=COVERAGE)


stratified = expand("results/stratified/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.extended.csv", tool=TOOL,  aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, bamtype=BAM, coverage=COVERAGE) + \
    expand("results/stratified/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.extended.csv", tool=TOOL,  aligner=ONT_ALIGNER, sample=ONT_SAMPLE, bamtype=BAM, coverage=COVERAGE)

snpeff = expand("results/SnpEff/{tool}_origin_{aligner}_{sample}_5.vcf.gz",  tool=TOOL,  aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, bamtype=BAM) + \
    expand("results/SnpEff/{tool}_origin_{aligner}_{sample}_5.vcf.gz",  tool=TOOL, aligner=ONT_ALIGNER, sample=ONT_SAMPLE, bamtype=BAM) 

ISOSEQ_SAMPLE = [s for s in ISOSEQ_SAMPLE if "HG002" in s or "HG005" in s]
MASSEQ_SAMPLE = [s for s in MASSEQ_SAMPLE if "HG002" in s or "HG005" in s]
ONT_SAMPLE    = [s for s in ONT_SAMPLE    if "HG002" in s or "HG005" in s]
TOOL = ["Clair3-RNA","DeepVariant", "GATK", "longcallR-nn"]
PACBIO_SAMPLE = ISOSEQ_SAMPLE + MASSEQ_SAMPLE

phase_bed = expand(
    "results/switch_error/{phaser}_{tool}_origin_minimap2_{sample}.bed",
    phaser=PHASER,
    tool=TOOL,
    sample=PACBIO_SAMPLE
) + expand(
    "results/switch_error/{phaser}_{tool}_origin_minimap2_{sample}.bed",
    phaser=PHASER,
    tool=TOOL,
    sample=ONT_SAMPLE
)

phase_tsv = expand(
    "results/switch_error/{phaser}_{tool}_origin_minimap2_{sample}.tsv",
    phaser=PHASER,
    tool=TOOL,
    sample=PACBIO_SAMPLE
) + expand(
    "results/switch_error/{phaser}_{tool}_origin_minimap2_{sample}.tsv",
    phaser=PHASER,
    tool=TOOL,
    sample=ONT_SAMPLE)

align = {
    "preprocess": preprocess,
    "origin": origin,
    "transformed": transformed,
}

eva = {
    #"qual": qual,
    "happy": happy,
    #"transformed": happy_transformed,
    "stratified": stratified,
    #"phase": phase_tsv,
    #"snpeff": snpeff
}

VAL = eva.keys()

plt = []
for val in VAL:
    x = glob_wildcards("code/plt-"+val+"_{x}.R").x
    plt += expand("plts/{val}-{plt}.pdf", val=val, plt=x)

# SETUP ------------------------------------------------------------------------------
rule all: 
    input: 
        #[x for x in align.values()], 
        [x for x in eva.values()],
        plt


############### Visualization #################
for val in VAL:
    rule:
        priority: 90
        input:  expand("code/plt-{val}_{{plt}}.R", val=val), x=eva[val]
        params: lambda wc, input: ";".join(input.x)
        output: expand("plts/{val}-{{plt}}.pdf", val=val)
        log:    expand("logs/plt_{val}-{{plt}}.Rout", val=val)
        shell:  '''
            {R} CMD BATCH --no-restore --no-save "--args\
            {params} {output[0]}" {input[0]} {log}'''
