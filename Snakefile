import json
import itertools
configfile: "config.yaml"
R = config["R"]

# WILDCARDS --------------------------------------------------------------------
ISOSEQ_SAMPLE = [
    "HG004-Baylor-IsoSeq",
    "HG005-Baylor-IsoSeq", 
    "HG002a-Baylor-IsoSeq"
   # "HG002-Baylor-IsoSeq"
]
MASSEQ_SAMPLE = ["HG004-MasSeq","HG002a-MasSeq", "HG005-MasSeq"]
PACBIO_SAMPLE = ISOSEQ_SAMPLE + MASSEQ_SAMPLE
ONT_SAMPLE = [
    "HG004-cDNAxR09", "HG004-cDNAxR10",  "HG004-dRNA004" ,"HG004-dRNA002",
    "HG005-cDNAxR09", "HG005-cDNAxR10",  "HG005-dRNA004" ,"HG005-dRNA002",
    "HG002a-cDNAxR10", "HG002a-dRNA002", "HG002a-cDNAxR09", "HG002-dRNA004"
    ] 

BAM = ["origin"] # transformation of bam
PACBIO_ALIGNER = ["minimap2", "pbmm2"] # pacbio aligner, minimap2-nojunc
ONT_ALIGNER = ["minimap2"] # ONT aligner, minimap2-nojunc
HG004 = [s for s in ONT_SAMPLE + PACBIO_SAMPLE if "HG004" in s]
COVERAGE = [1, 5, 10, 30, 50, 100] # coverage
TOOL = [ "Clair3-RNA",  "longcallR-nn", "DeepVariant", "longcallR",  "isoLASER", "GATK"] 
ISOSEQ_STATUS = ["raw", "preprocessed"]

# STRAT_BEDS = {}
# with open("data/regions/GRCh38-all-stratifications.tsv") as _fh:
#     for _line in _fh:
#         _name, _path = _line.rstrip("\n").split("\t")
#         STRAT_BEDS[_name] = "data/regions/" + _path

# CONTEXT_BINS = {
#     "GC content": [
#         "gc15_slop50", "gc15to20_slop50", "gc20to25_slop50", "gc25to30_slop50",
#         "gc30to55_slop50", "gc55to60_slop50", "gc60to65_slop50", "gc65to70_slop50",
#         "gc70to75_slop50", "gc75to80_slop50", "gc80to85_slop50", "gc85_slop50",
#     ],
#     "Homopolymer": [
#         "SimpleRepeat_homopolymer_4to6_slop5", "SimpleRepeat_homopolymer_7to11_slop5",
#         "SimpleRepeat_homopolymer_ge12_slop5", "SimpleRepeat_homopolymer_ge21_slop5",
#     ],
#     "Tandem Repeat": [
#         "AllTandemRepeats_le50bp_slop5", "AllTandemRepeats_51to200bp_slop5",
#         "AllTandemRepeats_201to10000bp_slop5", "AllTandemRepeats_ge10001bp_slop5",
#     ],
# }
# CONTEXT_SUBSETS = [s for bins in CONTEXT_BINS.values() for s in bins]

# RULES --------------------------------------------------------------------
include: "workflow/rules/align.smk"
include: "workflow/rules/variant_calling.smk"
include: "workflow/rules/variant_benchmark.smk"
include: "workflow/rules/phasing.smk"
include: "workflow/rules/phasing_benchmark.smk"


# RESULTS --------------------------------------------------------------------

#### Alignment 
preprocess = expand("data/preprocessed/{sample}.bam", sample=ISOSEQ_SAMPLE)

origin = expand("results/align/{status}/origin/{aligner}_{sample}.aligned.bam",
                aligner=PACBIO_ALIGNER, sample=ISOSEQ_SAMPLE, status=ISOSEQ_STATUS) + \
        expand("results/align/raw/origin/{aligner}_{sample}.aligned.bam",
                aligner=PACBIO_ALIGNER, sample=MASSEQ_SAMPLE) + \
        expand("results/align/raw/origin/{aligner}_{sample}.aligned.bam",
                aligner=ONT_ALIGNER, sample=ONT_SAMPLE) + \
        expand("results/align/raw/origin/{aligner}_{sample}.aligned.bam",
                  aligner="minimap2-nojunc", sample=HG004)

transformed = expand("results/align/{status}/transformed/minimap2_{sample}.aligned.bam", 
    sample=ISOSEQ_SAMPLE, status=ISOSEQ_STATUS) + \
    expand("results/align/raw/transformed/minimap2_{sample}.aligned.bam",
             sample=MASSEQ_SAMPLE) + \
    expand("results/align/raw/transformed/minimap2_{sample}.aligned.bam", sample=ONT_SAMPLE)

qual = expand("results/align/{status}/quality/{aligner}_{sample}.tsv", aligner=PACBIO_ALIGNER, sample=ISOSEQ_SAMPLE, status=ISOSEQ_STATUS) + \
    expand("results/align/raw/quality/{aligner}_{sample}.tsv", aligner=PACBIO_ALIGNER, sample=MASSEQ_SAMPLE) + \
    expand("results/align/raw/quality/{aligner}_{sample}.tsv", aligner=ONT_ALIGNER, sample=ONT_SAMPLE) 

bases = expand("results/beds/raw/summary/{sample}_{coverage}.tsv", sample = PACBIO_SAMPLE + ONT_SAMPLE, coverage = COVERAGE)

reads = expand("results/align/{status}/supplementary_stats/{aligner}_{sample}.tsv", aligner=PACBIO_ALIGNER, sample=ISOSEQ_SAMPLE, status=ISOSEQ_STATUS) + \
    expand("results/align/raw/supplementary_stats/{aligner}_{sample}.tsv", aligner=PACBIO_ALIGNER, sample=MASSEQ_SAMPLE) + \
    expand("results/align/raw/supplementary_stats/{aligner}_{sample}.tsv", aligner="minimap2-nojunc", sample=HG004) 

happy = expand("results/happy/{tool}_origin_{aligner}_{sample}_{coverage}.summary.csv", tool=TOOL, aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, coverage=COVERAGE) + \
   expand("results/happy/{tool}_origin_{aligner}_{sample}_{coverage}.summary.csv", tool=TOOL,  aligner=ONT_ALIGNER, sample=ONT_SAMPLE, coverage=COVERAGE) + \
    expand("results/happy/{tool}_origin_{aligner}_{sample}_{coverage}.summary.csv", tool=TOOL,  aligner=ONT_ALIGNER, sample=ONT_SAMPLE, coverage=COVERAGE) + \
    expand("results/happy/{tool}_origin_minimap2-nojunc_{sample}_{coverage}.summary.csv",
       tool=[t for t in TOOL if t != "GATK"], sample=HG004, coverage=COVERAGE)


happy_transformed = expand("results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.summary.csv", tool=TOOL, aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, coverage=COVERAGE, bamtype=BAM) + \
   expand("results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.summary.csv", tool=TOOL,  aligner=ONT_ALIGNER, sample=ONT_SAMPLE, coverage=COVERAGE, bamtype=BAM)


stratified = expand("results/stratified/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.extended.csv", tool=TOOL,  aligner="minimap2", sample=PACBIO_SAMPLE, bamtype=BAM, coverage=5) + \
    expand("results/stratified/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.extended.csv", tool=TOOL,  aligner="minimap2", sample=ONT_SAMPLE, bamtype=BAM, coverage=5)

# context_bases = expand("results/beds/context/{sample}_{coverage}_{subset}.tsv",
#     sample=HG004, coverage=5, subset=CONTEXT_SUBSETS)

snpeff = expand("results/SnpEff/{tool}_origin_{aligner}_{sample}_5.vcf.gz",  tool=TOOL,  aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, bamtype=BAM) + \
    expand("results/SnpEff/{tool}_origin_{aligner}_{sample}_5.vcf.gz",  tool=TOOL, aligner=ONT_ALIGNER, sample=ONT_SAMPLE, bamtype=BAM)

rna_editing = expand("results/rna_editing/{tool}_origin_minimap2_{sample}_{coverage}.tsv", tool=TOOL, sample=PACBIO_SAMPLE + ONT_SAMPLE, coverage=COVERAGE) + \
      expand("results/rna_editing/spectrum/{tool}_origin_minimap2_{sample}_{coverage}.tsv", tool=TOOL, sample=PACBIO_SAMPLE + ONT_SAMPLE, coverage=COVERAGE) #+ \

happy_chr20 = expand("results/happy_chr20/{tool}_origin_minimap2_{sample}_{coverage}.tsv",
    tool=TOOL, sample=["HG004-MasSeq", "HG005-MasSeq"], coverage=COVERAGE)

ISOSEQ_SAMPLE = [s for s in ISOSEQ_SAMPLE if "HG002" in s or "HG005" in s]
MASSEQ_SAMPLE = [s for s in MASSEQ_SAMPLE if "HG002" in s or "HG005" in s]
ONT_SAMPLE    = [s for s in ONT_SAMPLE    if "HG002" in s or "HG005" in s]
PHASER = ["HiPhase", "WhatsHap", "longphase", "HapCUT2"] 
TOOL = ["Clair3-RNA","DeepVariant", "longcallR-nn"]
PACBIO_SAMPLE = ISOSEQ_SAMPLE + MASSEQ_SAMPLE
SAMPLE = PACBIO_SAMPLE + ONT_SAMPLE


align = {
    #"preprocess": preprocess,
    "origin": origin,
    #"transformed": transformed,
}

TYPE = ["all", "SNP"]          
phase_tsv = (
    expand(
        "results/switch_error/summary/{phaser}_{tool}_origin_minimap2_{sample}_{type}.tsv",
        phaser=[p for p in PHASER if p != "HiPhase"],
        tool=TOOL,
        sample=SAMPLE,
        type=TYPE
    )
    + \
    expand(
        "results/switch_error/summary/HiPhase_{tool}_origin_minimap2_{sample}_{type}.tsv",
        tool=TOOL,
        sample=PACBIO_SAMPLE,
        type=TYPE
    ) + \
    expand(
        "results/switch_error/summary/longcallR_longcallR_origin_minimap2_{sample}_{type}.tsv",
        sample=SAMPLE,
        type=TYPE
    ) + \
    expand(
        "results/switch_error/summary/isoLASER_isoLASER_origin_minimap2_{sample}_{type}.tsv",
        sample=SAMPLE,
        type=TYPE
    )

)

phase_block_tsv = (
    expand(
        "results/phase_blocks/summary/{phaser}_{tool}_origin_minimap2_{sample}_{type}.tsv",
        phaser=[p for p in PHASER if p != "HiPhase"],
        tool=TOOL,
        sample=SAMPLE,
        type=TYPE
    )
    + \
    expand(
        "results/phase_blocks/summary/HiPhase_{tool}_origin_minimap2_{sample}_{type}.tsv",
        tool=TOOL,
        sample=PACBIO_SAMPLE,
        type=TYPE
    ) + \
    expand(
        "results/phase_blocks/summary/longcallR_longcallR_origin_minimap2_{sample}_{type}.tsv",
        sample=SAMPLE,
        type=TYPE
    ) + \
    expand(
        "results/phase_blocks/summary/isoLASER_isoLASER_origin_minimap2_{sample}_{type}.tsv",
        sample=SAMPLE,
        type=TYPE
    )

)

eva = {
    #"qual": qual,
    "bases": bases,
    "reads": reads,
    "happy": happy,
    "rna_editing": rna_editing,
    "happy-chr20": happy_chr20,
    #"transformed": happy_transformed,
    "stratified": stratified,
    #"phase": phase_tsv,
    #"phaseblocks": phase_block_tsv,
    "phasecombined": phase_tsv + phase_block_tsv,
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
        [x for x in align.values()], 
        [x for x in eva.values()],
        plt


############### Visualization #################
for val in VAL:
    rule:
        priority: 90
        conda:  "workflow/envs/r_env.yaml"
        input:  expand("code/plt-{val}_{{plt}}.R", val=val), x=eva[val]
        params: lambda wc, input: ";".join(input.x)
        output: expand("plts/{val}-{{plt}}.pdf", val=val)
        log:    expand("logs/plt_{val}-{{plt}}.Rout", val=val)
        shell: "R CMD BATCH --no-restore --no-save '--args {params} {output[0]}' {input[0]} {log[0]}"
