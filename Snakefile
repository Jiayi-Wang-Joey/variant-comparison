import json
import itertools
configfile: "config.yaml"
R = config["R"]

# WILDCARDS --------------------------------------------------------------------
ISOSEQ_SAMPLE = [#"WTC11-IsoSeq", "Jurkat-IsoSeq", "HG004-Baylor-IsoSeq", 
"HG004-Baylor-IsoSeq", "HG005-Baylor-IsoSeq"]
MASSEQ_SAMPLE = ["HG004-MasSeq","HG005-MasSeq"]
PACBIO_SAMPLE = ISOSEQ_SAMPLE + MASSEQ_SAMPLE
ONT_SAMPLE = ["HG004-dRNA002", "HG005-dRNA002", "HG004-cDNAxR09", 
    "HG005-cDNAxR09", "HG004-dRNA004", "HG005-dRNA004" ] 

BAM = ["origin", "transformed"] # transformation of bam
PACBIO_ALIGNER = ["minimap2", "pbmm2"] # pacbio aligner
ONT_ALIGNER = ["minimap2"] # ONT aligner
COVERAGE = [5, 10, 30, 50, 100] # coverage
TOOL = ["Clair3-RNA", "longcallR", "DeepVariant"] # variant caller
PHASER = ["WhatsHap", "longphase"]

ISOSEQ_STATUS = ["raw", "preprocessed"]

# RESULTS ------------------------------------------------------------------------------

# alignment
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

# data quality 

qual = expand("results/align/{status}/quality/{aligner}_{sample}.tsv", aligner=PACBIO_ALIGNER, sample=ISOSEQ_SAMPLE, status=ISOSEQ_STATUS) + \
    expand("results/align/raw/quality/{aligner}_{sample}.tsv", aligner=PACBIO_ALIGNER, sample=MASSEQ_SAMPLE) + \
    expand("results/align/raw/quality/{aligner}_{sample}.tsv", aligner=ONT_ALIGNER, sample=ONT_SAMPLE) 

cov = expand("results/align/{status}/coverage/{aligner}_{sample}.tsv", aligner=PACBIO_ALIGNER, sample=ISOSEQ_SAMPLE, status=ISOSEQ_STATUS) + \
    expand("results/align/raw/coverage/{aligner}_{sample}.tsv", aligner=PACBIO_ALIGNER, sample=MASSEQ_SAMPLE) + \
    expand("results/align/raw/coverage/{aligner}_{sample}.tsv", aligner=ONT_ALIGNER, sample=ONT_SAMPLE) 

# variant calling
clair3_rna = expand("results/variant/Clair3-RNA/{bamtype}_{aligner}_{sample}.vcf.gz", 
    aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, bamtype=BAM) + \
    expand("results/variant/Clair3-RNA/{bamtype}_{aligner}_{sample}.vcf.gz", 
    aligner=ONT_ALIGNER, sample=ONT_SAMPLE, bamtype=BAM)

deep_varaint = expand("results/variant/DeepVariant/{bamtype}_{aligner}_{sample}.vcf.gz", 
    aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, bamtype=BAM) + \
    expand("results/variant/DeepVariant/{bamtype}_{aligner}_{sample}.vcf.gz", 
    aligner=ONT_ALIGNER, sample=ONT_SAMPLE, bamtype=BAM)

longcallR = directory(expand("results/variant/longcallR/{bamtype}_{aligner}_{sample}", 
    aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, bamtype=BAM)) + \
    directory(expand("results/variant/longcallR/{bamtype}_{aligner}_{sample}", 
    aligner=ONT_ALIGNER, sample=ONT_SAMPLE, bamtype=BAM))


# benchmarking
bed = expand("results/beds/{status}/{sample}_{aligner}_{coverage}.bed", aligner=PACBIO_ALIGNER, sample=ISOSEQ_SAMPLE, coverage=COVERAGE, status=ISOSEQ_STATUS) + \
    expand("results/beds/raw/{sample}_{aligner}_{coverage}.bed", aligner=ONT_ALIGNER, sample=ONT_SAMPLE, coverage=COVERAGE)

happy = expand("results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.summary.csv", tool=TOOL, aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, bamtype=BAM, coverage=COVERAGE) + \
   expand("results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.summary.csv", tool=TOOL,  aligner=ONT_ALIGNER, sample=ONT_SAMPLE, bamtype=BAM, coverage=COVERAGE)

roc = expand("results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.roc.csv.gz",  tool=TOOL, aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, bamtype=BAM,  coverage=COVERAGE) + \
    expand("results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.roc.csv.gz",  tool=TOOL, aligner=ONT_ALIGNER, sample=ONT_SAMPLE, bamtype=BAM,  coverage=COVERAGE)


stratified = expand("results/stratified/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.extended.csv", tool=TOOL,  aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, bamtype=BAM, coverage=COVERAGE) + \
    expand("results/stratified/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.extended.csv", tool=TOOL,  aligner=ONT_ALIGNER, sample=ONT_SAMPLE, bamtype=BAM, coverage=COVERAGE)

## phasing
ISOSEQ_SAMPLE = [s for s in ISOSEQ_SAMPLE if "HG002" in s or "HG005" in s]
MASSEQ_SAMPLE = [s for s in MASSEQ_SAMPLE if "HG002" in s or "HG005" in s]
ONT_SAMPLE    = [s for s in ONT_SAMPLE    if "HG002" in s or "HG005" in s]

PACBIO_SAMPLE = ISOSEQ_SAMPLE + MASSEQ_SAMPLE

# phase = expand("results/phase/{phaser}/{tool}_{bamtype}_{aligner}_{sample}.vcf.gz", phaser = PHASER, tool = TOOL, bamtype=BAM, aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE) + \
# expand("results/phase/{phaser}/{tool}_{bamtype}_{aligner}_{sample}.vcf.gz", phaser = PHASER, tool = TOOL, bamtype=BAM, aligner=ONT_ALIGNER, sample=ONT_SAMPLE) 
phase_bed = expand("results/switch_error/{phaser}_{tool}_{bamtype}_{aligner}_{sample}.bed", phaser = PHASER, tool = TOOL, bamtype=BAM, aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE) + \
expand("results/switch_error/{phaser}_{tool}_{bamtype}_{aligner}_{sample}.bed", phaser = PHASER, tool = TOOL, bamtype=BAM, aligner=ONT_ALIGNER, sample=ONT_SAMPLE) 
phase_tsv =  expand("results/switch_error/{phaser}_{tool}_{bamtype}_{aligner}_{sample}.tsv", phaser = PHASER, tool = TOOL, bamtype=BAM, aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE) + \
expand("results/switch_error/{phaser}_{tool}_{bamtype}_{aligner}_{sample}.tsv", phaser = PHASER, tool = TOOL, bamtype=BAM, aligner=ONT_ALIGNER, sample=ONT_SAMPLE) 

res = {
    # "origin": origin,
    # "transformed": transformed,
    # "qual": qual,
    # "cov": cov,
    # "clair3-rna": clair3_rna,
    # "deep_variant": deep_varaint,
    # "longcallR": longcallR,
    # "bed": bed,
    "happy": happy,
    # "roc": roc,
    # "stratified": stratified,
    "phase": phase_tsv
}

# visualization
VAL = res.keys()

plt = []
for val in VAL:
    x = glob_wildcards("code/plt-"+val+"_{x}.R").x
    plt += expand("plts/{val}-{plt}.pdf", val=val, plt=x)

# ------------------------------------------------------------------------------
rule all: 
    input: 
        [x for x in res.values()], plt

############################ Alignment ##################################

#rule get_data:

rule preprocess: 
    priority: 100
    input: 
        bam="data/raw/{sample}.bam",
        primer="data/primer/primer_isoseq.fasta"
    output: "data/preprocessed/{sample}.bam"
    conda: "envs/preprocessing.yaml"
    log:
        stderr="logs/preprocess/{sample}.stderr"
    wildcard_constraints:
        sample=".*IsoSeq.*"
    shell: 
        """
        lima {input.bam} {input.primer} data/raw/{wildcards.sample}.fl.bam --isoseq --peek-guess 2> {log}    
        isoseq refine data/raw/{wildcards.sample}.fl.NEB_5p--NEB_Clontech_3p.bam {input.primer} {output} 2> {log}        
        """       

def extract_tech(wildcards):
    return wildcards.sample.split("-")[-1] 

rule bam2fq:
    priority: 99
    input:
        "data/{status}/{sample}.bam"
    output:
        "results/align/{status}/bam2fq/{sample}.fastq.gz"
    log:
        stderr="logs/align/{status}/bam2fq/{sample}.stderr"
    threads: config["convert_threads"]
    params:
        threads=config["convert_threads"],
    wildcard_constraints:
        status="raw|preprocessed" if "{sample}" in ISOSEQ_SAMPLE else "raw" 
    shell:
        """
        samtools bam2fq -@ {params.threads} {input} | gzip > {output}
        """

def minimap2_preset(wildcards):
    if "IsoSeq" in wildcards.sample:
        return "-ax splice:hq -uf"
    elif "MasSeq" in wildcards.sample:
        return "-ax splice:hq -uf"
    elif "dRNA" in wildcards.sample:
        return "-ax splice -uf -k14"
    elif "cDNA" in wildcards.sample:
        return "-ax splice"
    else:
        raise ValueError(f"Cannot determine minimap2 preset from sample name: {wildcards.sample}")


rule minimap2_align:
    priority: 98
    input:
        reads = rules.bam2fq.output,
        transcriptome = config["transcriptome_bed"],
        genome = config["reference_genome"],
    output:
        "results/align/{status}/origin/minimap2_{sample}.aligned.bam"
    params:
        align_map_bam_threads = config["align_map_bam_threads"],
        align_sort_bam_threads = config["align_sort_bam_threads"],
        align_sort_bam_memory_gb = config["align_sort_bam_memory_gb"],
        preset = minimap2_preset
    conda:
        "envs/align.yaml"
    log:
        stdout = "logs/minimap2/{status}/{sample}.out",
        stderr = "logs/minimap2/{status}/{sample}.err"
    wildcard_constraints:
        status="raw|preprocessed" if "{sample}" in ISOSEQ_SAMPLE else "raw" 
    shell:
        """
        minimap2 {params.preset} --junc-bed {input.transcriptome} \
            -t {params.align_map_bam_threads} \
            {input.genome} {input.reads} | samtools sort \
            -@ {params.align_sort_bam_threads} \
            -m{params.align_sort_bam_memory_gb}g \
            -o {output} > {log.stdout} 2> {log.stderr}
        """


rule pbmm2_align:
    priority: 98
    input:
        reads = "results/align/{status}/bam2fq/{sample}.fastq.gz",   
        genome = config["reference_genome"]
    output:
        bam = "results/align/{status}/origin/pbmm2_{sample}.aligned.bam"
    params:
        preset = "CCS",  
        threads = config["align_map_bam_threads"]
    conda:
        "envs/align.yaml"
    log:
        stdout = "logs/pbmm2/{status}/{sample}.out",
        stderr = "logs/pbmm2/{status}/{sample}.err"
    wildcard_constraints:
        status="raw|preprocessed" if "{sample}" in ISOSEQ_SAMPLE else "raw" 
    shell:
        """
        pbmm2 align \
            --preset {params.preset} \
            --sort \
            --num-threads {params.threads} \
            {input.genome} \
            {input.reads} \
            {output.bam} > {log.stdout} 2> {log.stderr}
        """


rule graphmap2_align:
    input:
        reads = "results/align/genome/bam2fq/{sample}.fastq.gz",   
        genome = config["reference_genome"]
    output:
        "results/align/genome/origin/GraphMap2_{sample}.aligned.bam"
    params:
        threads = config["align_map_bam_threads"],
        sort_threads = config["align_sort_bam_threads"],
        sort_memory = config["align_sort_bam_memory_gb"]
    log:
        stdout = "logs/graphmap2/{sample}.out",
        stderr = "logs/graphmap2/{sample}.err"
    shell:
        """
        ~/tools/graphmap2/bin/Linux-x64/graphmap2 align \
            -x rnaseq \
            -r {input.genome} \
            -d {input.reads} \
            -t {params.threads} | \
        samtools sort \
            -@ {params.sort_threads} \
            -m {params.sort_memory}g \
            -o {output} > {log.stdout} 2> {log.stderr}
        """


def ultra_preset(wildcards):
    if "IsoSeq" in wildcards.sample or "MasSeq" in wildcards.sample:
        return "--isoseq"
    elif "dRNA" in wildcards.sample or "cDNA" in wildcards.sample:
        return "--ont"
    else:
        raise ValueError(f"Cannot determine minimap2 preset from sample name: {wildcards.sample}")


rule uLTRA_align:
    input:
        reads = "results/align/genome/bam2fq/{sample}.fastq.gz",   
        genome = config["reference_genome"]
    output:
        "results/align/genome/origin/uLTRA_{sample}.aligned.bam"
    params:
        threads = config["align_map_bam_threads"],
        sort_threads = config["align_sort_bam_threads"],
        sort_memory = config["align_sort_bam_memory_gb"],
        preset = ultra_preset
    log:
        stdout = "logs/uLTRA/{sample}.out",
        stderr = "logs/uLTRA/{sample}.err"
    shell:
        """
        mkdir -p results/align/genome/origin/ultra/

        uLTRA align \
            {input.genome} \
            {input.reads} \
            results/align/genome/origin/ultra/ \
            {params.preset} \
            --t {params.threads} \
            --prefix {wildcards.sample} > {log.stdout} 2> {log.stderr}

        samtools sort \
            -@ {params.sort_threads} \
            -m {params.sort_memory}G \
            -o {output} \
            results/align/genome/origin/ultra/{wildcards.sample}.sam

        rm results/align/genome/origin/ultra/{wildcards.sample}.sam
        """


############################ Data quality ####################################

rule bam_quality:
    priority: 80
    input:
        "results/align/{status}/origin/{aligner}_{sample}.aligned.bam"
    output:
        stat = "results/align/{status}/quality/{aligner}_{sample}.tsv",
        cov = "results/align/{status}/coverage/{aligner}_{sample}.tsv"
    log:
        stdout = "logs/quality/{status}_{aligner}_{sample}.log",
        stderr = "logs/quality/{status}_{aligner}_{sample}.err"
    shell:
        """
        echo -e "sample\\tstatus\\taligner\\tmetric\\tvalue\\tannotation" > {output.stat}
        echo -e "sample\\tstatus\\taligner\\trange\\tdepth\\tcount" > {output.cov}

        samtools stats {input} 2>> {log.stderr} | tee >(grep ^SN | \
            awk -v sample="{wildcards.sample}" -v status="{wildcards.status}" -v aligner="{wildcards.aligner}" 'BEGIN {{FS="\\t"; OFS="\\t"}} \
            {{metric = $2; value = $3; annotation = ($0 ~ /#/) ? substr($0, index($0, "#")) : ""; print sample, status, aligner, metric, value, annotation}}' >> {output.stat}) \
        | grep ^COV | cut -f 2- | \
            awk -v s="{wildcards.sample}" -v t="{wildcards.status}" -v a="{wildcards.aligner}" 'BEGIN {{OFS="\\t"}} {{print s,t,a,$1,$2,$3}}' >> {output.cov}
        """


############################ Transformation ##################################
rule split_ncigar:
    priority: 97
    input:
        bam = "results/align/{status}/origin/{aligner}_{sample}.aligned.bam",
        ref = config["reference_genome"],
        img = "/home/jiayiwang/tools/gatk-latest.simg"  
    output: 
        "results/align/{status}/splitNC/{aligner}_{sample}.aligned.bam"
    log:
        "logs/split_ncigar_reads_{status}_{aligner}_{sample}.log"
    params:
        threads = 10
    wildcard_constraints:
        status="raw|preprocessed" if "{sample}" in ISOSEQ_SAMPLE else "raw" 
    shell:
        """
        singularity exec --bind /home/jiayiwang/miniconda3/envs/snakemake/bin/python:/usr/bin/python {input.img} /gatk/gatk --java-options \
        "-Xmx16G -XX:+UseParallelGC -XX:ParallelGCThreads={params.threads} -Djava.util.concurrent.ForkJoinPool.common.parallelism={params.threads} -Djava.io.tmpdir=/home/jiayiwang/tmp1" SplitNCigarReads \
            -R {input.ref} \
            -I {input.bam} \
            -O {output} \
        > {log} 2>&1
        """ 

rule flag_correction:
    priority: 96
    input:
        bam = "results/align/{status}/origin/{aligner}_{sample}.aligned.bam",
        sncr_bam = rules.split_ncigar.output
    output: 
        "results/align/{status}/transformed/{aligner}_{sample}.aligned.bam"
    log:
        "logs/flag_correction_{status}_{aligner}_{sample}.log"
    params:
        threads = 5
    wildcard_constraints:
        status="raw|preprocessed" if "{sample}" in ISOSEQ_SAMPLE else "raw" 
    shell:
        """
        Rscript /home/jiayiwang/tools/lrRNAseqVariantCalling/tools/flagCorrection.r \
            {input.bam} \
            {input.sncr_bam} \
            {output} \
            {params.threads} \
        > {log} 2>&1
        rm {input.sncr_bam}
        """

############################ Variant Calling ##################################

def get_clair3_rna_platform(wildcards):
    sample = wildcards.sample.lower()
    aligner = wildcards.aligner.lower()

    if "masseq" in sample:
        if aligner in ["pbmm2", "minimap2"]:
            return f"hifi_mas_{aligner}"
        else:
            return "hifi_mas"
    elif "isoseq" in sample:
        if aligner in ["pbmm2", "minimap2"]:
            return f"hifi_sequel2_{aligner}"
        else:
            return "hifi_sequel2"
    elif "drna004" in sample:
        return "ont_dorado_drna004"
    elif "r10" in sample and "cdna" in sample:
        return "ont_r10_dorado_cdna"
    elif "drna002" in sample:
        return "ont_guppy_drna002"
    elif "r09" in sample and "cdna" in sample:
        return "ont_r9_guppy_cdna"
    else:
        raise ValueError(f"Cannot determine platform for sample '{wildcards.sample}' and aligner '{wildcards.aligner}'")

rule run_clair3_rna:
    priority: 95
    input:
        bam="results/align/raw/{bamtype}/{aligner}_{sample}.aligned.bam",
        ref=config["reference_genome"],
        img="/home/jiayiwang/tools/clair3-rna-latest.simg"
    output:
        dir=directory("/home/jiayiwang/variant-comparison/results/variant/Clair3-RNA/{bamtype}_{aligner}_{sample}"),
        file="results/variant/Clair3-RNA/{bamtype}_{aligner}_{sample}.vcf.gz"
    log: 
        "logs/run_clair3_rna_{bamtype},{aligner},{sample}.log"
    params:
        threads = 10,
        model = get_clair3_rna_platform
    shell:
        """
        if [ ! -f {input[0]}.bai ]; then
            samtools index {input[0]};
        fi
        singularity exec -B /home/jiayiwang/variant-comparison/{input[0]},{input[1]} {input[2]} \
            /bin/bash -c "source /opt/conda/bin/activate /opt/conda/envs/clair3_rna && \
            /opt/bin/run_clair3_rna \
            --bam_fn /home/jiayiwang/variant-comparison/{input[0]} \
            --ref_fn {input[1]} \
            --threads {params.threads} \
            --platform {params.model} \
            --tag_variant_using_readiportal \
            --remove_intermediate_dir \
            --output_dir {output[0]} \
            --conda_prefix /opt/conda/envs/clair3_rna" \
            > {log} 2>&1
        mv {output[0]}/output.vcf.gz {output[1]}
        """

# rule run_clair3:
#     priority: 95
#     input:
#         bam="/home/jiayiwang/variant-comparison/results/align/genome/{bamtype}/{aligner}_{sample}.aligned.bam",
#         ref=config["reference_genome"],
#         img="/home/jiayiwang/tools/clair3-latest.simg"
#     output:
#         dir=directory("/home/jiayiwang/variant-comparison/results/variant/Clair3/{bamtype}_{aligner}_{sample}"),
#         file="/home/jiayiwang/variant-comparison/results/variant/Clair3/{bamtype}_{aligner}_{sample}.vcf.gz"
#     log: 
#         "logs/run_clair3/{bamtype},{aligner},{sample}.log"
#     params:
#         threads = 10,
#         model = get_clair3_platform
#     shell:
#         """
#         if [ ! -f {input[0]}.bai ]; then
#             samtools index {input[0]};
#         fi
#         singularity exec -B {input[0]},{input[1]} 
#             {input.img} \
#             /opt/bin/run_clair3.sh \
#             --bam_fn={input[0]} \
#             --ref_fn={input[1]} \
#             --threads={params.threads} \
#             --platform={params.model} \
#             --output_dir {output[0]} \
#             --model_path="/opt/models/{params.model}" \
#             > {log} 2>&1
#         mv {output[0]}/output.vcf.gz {output[1]}
#         """
        

# rule run_deep_variant:
#     priority: 95
#     input: 
#         bam=lambda wildcards: 
#             f"results/align/preprocessed/{wildcards.bamtype}/{wildcards.aligner}_{wildcards.sample}.aligned.bam"
#             if wildcards.sample in ISOSEQ_SAMPLE
#             else f"results/align/raw/{wildcards.bamtype}/{wildcards.aligner}_{wildcards.sample}.aligned.bam",
#         ref=config["reference_genome"],
#         img="/home/jiayiwang/tools/deepvariant_1.9.0.sif"
#     output:
#         "results/variant/DeepVariant/{bamtype}_{aligner}_{sample}.vcf.gz",
#     log:
#         "logs/run_deep_variant_{bamtype},{aligner},{sample}.log"
#     threads: 10
#     params: 
#         threads = 10,
#         model = lambda wildcards: "MASSEQ" if "IsoSeq" in wildcards.sample or "MasSeq" in wildcards.sample else "ONT_R104",
#         extra_args = lambda wildcards: "" if "IsoSeq" in wildcards.sample or "MasSeq" in wildcards.sample else '--make_examples_extra_args="split_skip_reads=true"'
#     shell:
#         """
#         if [ ! -f {input.bam}.bai ]; then
#         samtools index {input.bam};
#         fi
#         singularity exec --bind /usr/lib/locale/ {input.img} /opt/deepvariant/bin/run_deepvariant\
#             --model_type {params.model} \
#             --ref {input.ref} \
#             --reads {input.bam} \
#             --output_vcf {output} \
#             --num_shards {threads} \
#             --intermediate_results_dir /home/jiayiwang/tmp/{wildcards.sample} \
#         > {log} 2>&1
#         """

def get_longcallR_platform(wildcards):
    sample = wildcards.sample.lower()
    for key, platform in {
        "masseq": "hifi-masseq",
        "isoseq": "hifi-isoseq",
        "cdna": "ont-cdna",
        "drna": "ont-drna"
    }.items():
        if key in sample:
            return platform
    raise ValueError(f"Unrecognized platform for sample '{wildcards.sample}'")


rule run_longcallR:
    priority: 95
    input: 
        "results/align/raw/{bamtype}/{aligner}_{sample}.aligned.bam"
    output:
        dir=directory("results/variant/longcallR/{bamtype}_{aligner}_{sample}"),
        file="results/variant/longcallR/{bamtype}_{aligner}_{sample}.vcf.gz"
    log:
        "logs/run_longcallR_{bamtype},{aligner},{sample}.log"
    params: 
        threads = 10,
        platform = get_longcallR_platform,
        ref = config["reference_genome"]
    shell:
        """
        if [ ! -f {input}.bai ] || [ {input} -nt {input}.bai ]; then
            samtools index {input};
        fi

        mkdir -p {output.dir}

        ~/tools/longcallR/target/release/longcallR \
            --bam-path {input} \
            --ref-path {params.ref} \
            --output {output.dir}/output \
            --preset {params.platform} \
            -t {params.threads} \
            --no-bam-output \
            --min-depth 0 \
            --max-depth 1000000 \
            > {log} 2>&1

        bgzip {output.dir}/output.vcf
        mv {output.dir}/output.vcf.gz {output.file}
        """

# rule run_longcallR-nn:
#     input: 
#         "results/align/raw/{bamtype}/{aligner}_{sample}.aligned.bam"
    
############################ Benchmarking ##################################
rule generate_bed:
    priority: 93
    input: 
        bam="results/align/{status}/origin/{aligner}_{sample}.aligned.bam",
        gtf="/home/jiayiwang/pacbio/data/genome/gencode.v47.chr_patch_hapl_scaff.annotation.bed"
    output:
        "results/beds/{status}/{sample}_{aligner}_{coverage}.bed"
    params:
        threads = 8,
        img = "/home/jiayiwang/tools/mosdepth-0.3.3--h37c5b7d_2.simg",
        min_cov=lambda w: w.coverage,
        confidence = lambda wildcards: f"data/truth/beds/{wildcards.sample.split('-')[0]}.bed",
    shell:
        """
        if [ ! -f {input[0]}.bai ]; then
            samtools index {input[0]};
        fi
        mkdir -p results/beds/{wildcards.status}/{wildcards.sample}
        singularity exec {params.img} mosdepth \
            --threads {params.threads} \
            --fasta config["reference_genome"] \
            results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage} \
            {input.bam}

        gzip -fdc results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.per-base.bed.gz \
            | awk -v min_cov={params.min_cov} '$4 >= min_cov' \
            | bedtools merge -d 1 -c 4 -o mean -i - \
            > results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.coverage.bed

        bedtools intersect -a results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.coverage.bed \
            -b {input.gtf} > results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.tmp.bed
        

        if [ -f {params.confidence} ]; then
            bedtools intersect -a results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.tmp.bed \
                -b {params.confidence} > {output}
            rm results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.tmp.bed
        else
            mv results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.tmp.bed {output}
        fi
        rm -rf results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}
        """

def get_bed(wildcards):
    sample = wildcards.sample
    tool = wildcards.tool
    aligner = wildcards.aligner
    coverage = wildcards.coverage

    if "IsoSeq" in sample and tool == "DeepVariant":
        status = "preprocessed"
    else:
        status = "raw"

    return f"results/beds/{status}/{sample}_{aligner}_{coverage}.bed"


rule happy_benchmark:
    priority: 92
    input:
        vcf="results/variant/{tool}/{bamtype}_{aligner}_{sample}.vcf.gz",
        bed=get_bed,
        confidence = lambda wildcards: f"data/truth/beds/{wildcards.sample.split('-')[0]}.bed"
    output:
        smy="results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.summary.csv",
        roc="results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.roc.csv.gz",
        vcf="results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.vcf.gz"
    params:
        truth = lambda wildcards: f"data/truth/{wildcards.sample.split('-')[0]}/truth.vcf.gz",
        ref = config["reference_genome"],
        outdir = "results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}",
        threads = 10 
    conda:
        "envs/happy.yaml"
    log:
        "logs/happy_{tool}_{bamtype}_{aligner}_{sample}_{coverage}.log"
    shell:
        """
            export HGREF={params.ref}
        
            if ! bcftools index -t -f "{input.vcf}"; then
                echo "VCF is not sorted. Sorting and reindexing..."
                cp {input.vcf} {input.vcf}.unsorted
                bcftools sort -Oz -o {input.vcf} {input.vcf}.unsorted
                bcftools index -t -f {input.vcf}
                rm {input.vcf}.unsorted
            fi
            singularity exec \
                --bind /home/jiayiwang/variant-comparison:/home/jiayiwang/variant-comparison \
                /home/jiayiwang/tools/hap.py.simg /opt/hap.py/bin/hap.py \
                {params.truth} \
                {input.vcf} \
                -r {params.ref} \
                -f {input.confidence} \
                -T {input.bed} \
                -o {params.outdir} \
                --engine=vcfeval \
                --threads={params.threads}
            
            awk -v tool="{wildcards.tool}" \
                -v bamtype="{wildcards.bamtype}" \
                -v aligner="{wildcards.aligner}" \
                -v sample="{wildcards.sample}" \
                -v coverage="{wildcards.coverage}" \
            'BEGIN {{FS=OFS=","}} 
                NR==1 {{print $0, "tool", "bamtype", "aligner", "sample", "coverage"}} 
                NR>1 {{print $0, tool, bamtype, aligner, sample, coverage}}' "{output.smy}" > "{output.smy}.tmp" && \
                mv "{output.smy}.tmp" "{output.smy}"

            zcat results/happy/{wildcards.tool}_{wildcards.bamtype}_{wildcards.aligner}_{wildcards.sample}_{wildcards.coverage}.roc.all.csv.gz | \
            awk -v tool="{wildcards.tool}" \
                -v bamtype="{wildcards.bamtype}" \
                -v aligner="{wildcards.aligner}" \
                -v sample="{wildcards.sample}" \
                -v coverage="{wildcards.coverage}" \
            'BEGIN {{FS=OFS=","}} 
                NR==1 {{print $0, "tool", "bamtype", "aligner", "sample", "coverage"}} 
                NR>1 {{print $0, tool, bamtype, aligner, sample, coverage}}'| gzip > {output.roc}
        """




rule stratified_benchmark:
    input:
        vcf=rules.happy_benchmark.output.vcf,
        stratification="/home/jiayiwang/variant-comparison/data/regions/v3.0-GRCh38-stratifications-all-except-genome-specific-stratifications.tsv"
    output:
        "results/stratified/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.extended.csv",
    params:
        ref=config["reference_genome"],
        outdir="results/stratified/{tool}_{bamtype}_{aligner}_{sample}_{coverage}",
        threads=5
    log:
        "logs/stratified_{tool}_{bamtype}_{aligner}_{sample}_{coverage}.log"
    shell:
        """
        singularity exec \
            --bind /home/jiayiwang/variant-comparison:/home/jiayiwang/variant-comparison \
            /home/jiayiwang/tools/hap.py.simg /opt/hap.py/bin/qfy.py \
            {input.vcf} \
            -o {params.outdir} \
            -t ga4gh \
            -r {params.ref} \
            --threads {params.threads} \
            --stratification {input.stratification}

        awk -v tool="{wildcards.tool}" -v bamtype="{wildcards.bamtype}" -v aligner="{wildcards.aligner}" -v sample="{wildcards.sample}" -v coverage="{wildcards.coverage}" '
        BEGIN {{FS=OFS=","}} 
        NR==1 {{print $0, "tool", "bamtype", "aligner", "sample", "coverage"}} 
        NR>1 {{print $0, tool, bamtype, aligner, sample, coverage}}' "{output}" > "{output}.tmp" && mv "{output}.tmp" "{output}"
        """

############### Phasing #####################

def get_bam(wildcards):
    sample = wildcards.sample
    tool = wildcards.tool
    aligner = wildcards.aligner

    if "IsoSeq" in sample and tool == "DeepVariant":
        status = "preprocessed"
    else:
        status = "raw"

    return f"results/align/{status}/origin/{aligner}_{sample}.aligned.bam"

rule WhatsHap:
    input:
        vcf="results/variant/{tool}/{bamtype}_{aligner}_{sample}.vcf.gz",
        bam=get_bam,
        ref=config["reference_genome"]
    output:
        "results/phase/WhatsHap/{tool}_{bamtype}_{aligner}_{sample}.vcf.gz"
    log:
        "logs/phase_WhatsHap_{tool}_{bamtype}_{aligner}_{sample}.log"
    shell:
        """
        whatshap phase --ignore-read-groups -o {output} --reference={input.ref} {input.vcf} {input.bam} > {log} 2>&1
        """

def get_platform_longphase(wildcards):
    sample = wildcards.sample
    if sample in PACBIO_SAMPLE:
        return "--pb"
    elif sample in ONT_SAMPLE:
        return "--ont"
    else:
        raise ValueError(f"Unknown platform for sample: {sample}")
    
rule longphase:
    input:
        vcf="results/variant/{tool}/{bamtype}_{aligner}_{sample}.vcf.gz",
        bam=get_bam,
        ref=config["reference_genome"]
    output:
        "results/phase/longphase/{tool}_{bamtype}_{aligner}_{sample}.vcf.gz"
    log:
        "logs/phase_longphase_{tool}_{bamtype}_{aligner}_{sample}.log"
    threads: 8
    params:
        platform=lambda wildcards: get_platform_longphase(wildcards)
    shell:
        """
        /home/jiayiwang/tools/longphase_linux-x64 phase \
            -s {input.vcf} \
            -b {input.bam} \
            -r {input.ref} \
            -t {threads} \
            -o results/phase/longphase/{wildcards.tool}_{wildcards.bamtype}_{wildcards.aligner}_{wildcards.sample} \
            --indels \
            {params.platform} \
            > {log} 2>&1

        bgzip results/phase/longphase/{wildcards.tool}_{wildcards.bamtype}_{wildcards.aligner}_{wildcards.sample}.vcf
        """

rule switch_error:
    input:
        vcf = "results/phase/{phaser}/{tool}_{bamtype}_{aligner}_{sample}.vcf.gz"
    output:
        bed = "results/switch_error/{phaser}_{tool}_{bamtype}_{aligner}_{sample}.bed",
        tsv = "results/switch_error/{phaser}_{tool}_{bamtype}_{aligner}_{sample}.tsv"
    log:
        "logs/phase_switch_error_{phaser}_{tool}_{bamtype}_{aligner}_{sample}.log"
    params:
        truth = lambda wildcards: f"data/truth/{wildcards.sample.split('-')[0]}/phase.vcf.gz"
    shell:
        """
        # Index truth VCF if missing
        if [ ! -f {params.truth}.tbi ]; then
            echo "Index file for truth VCF missing, creating with tabix..."
            tabix -p vcf {params.truth}
        fi

        # Index input VCF if missing
        if [ ! -f {input.vcf}.tbi ]; then
            echo "Index file for input VCF missing, creating with tabix..."
            tabix -p vcf {input.vcf}
        fi

        whatshap compare --ignore-sample-name \
            --switch-error-bed {output.bed} \
            --tsv-pairwise {output.tsv}.tmp \
            {params.truth} {input.vcf} > {log} 2>&1

        awk -v sample="{wildcards.sample}" \
            -v bamtype="{wildcards.bamtype}" \
            -v caller="{wildcards.tool}" \
            -v phaser="{wildcards.phaser}" \
            -v aligner="{wildcards.aligner}" \
            'BEGIN {{FS=OFS="\\t"}} NR==1 {{print "sample","bamtype","caller","phaser","aligner",$0}} NR>1 {{print sample,bamtype,caller,phaser,aligner,$0}}' \
            {output.tsv}.tmp > {output.tsv}

        rm {output.tsv}.tmp
        """




############### Visualization #################
for val in VAL:
    rule:
        priority: 90
        input:  expand("code/plt-{val}_{{plt}}.R", val=val), x=res[val]
        params: lambda wc, input: ";".join(input.x)
        output: expand("plts/{val}-{{plt}}.pdf", val=val)
        log:    expand("logs/plt_{val}-{{plt}}.Rout", val=val)
        shell:  '''
            {R} CMD BATCH --no-restore --no-save "--args\
            {params} {output[0]}" {input[0]} {log}'''
