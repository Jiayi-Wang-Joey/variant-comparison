import json
import itertools
configfile: "config.yaml"
R = config["R"]

# WILDCARDS --------------------------------------------------------------------
PACBIO_SAMPLE = ["WTC11-IsoSeq", "Jurkat-IsoSeq", "HG004-Baylor-IsoSeq", 
"HG005-Baylor-IsoSeq", "HG004-MasSeq"]
ONT_SAMPLE = ["HG005-dRNA002", "HG004-dRNA002"] # "HG004-dRNA002"
BAM = ["origin", "transformed"]
PACBIO_ALIGNER = ["minimap2", "pbmm2"]
ONT_ALIGNER = ["minimap2"]
# COVERAGE = ["1-5", "5-20", "20-50", "50-100", "gt100", 
# "gt5", "gt10", "gt30", "gt50"]
COVERAGE = ["gt5", "gt10", "gt30", "gt50", "gt100"]
TOOL = ["Clair3-RNA", "DeepVariant", "longcallR"]


# ------------------------------------------------------------------------------
# alignment
origin = expand("results/align/genome/origin/{aligner}_{sample}.aligned.bam", aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE) + \
    expand("results/align/genome/origin/{aligner}_{sample}.aligned.bam", aligner=ONT_ALIGNER, sample=ONT_SAMPLE)

transformed = expand("results/align/genome/transformed/{aligner}_{sample}.aligned.bam", aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE) + \
              expand("results/align/genome/transformed/{aligner}_{sample}.aligned.bam", aligner=ONT_ALIGNER, sample=ONT_SAMPLE)


# variant calling
clair3_rna = expand("/home/jiayiwang/variant-comparison/results/variant/Clair3-RNA/{bamtype}_{aligner}_{sample}.vcf.gz", 
    aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, bamtype=BAM) + \
    expand("/home/jiayiwang/variant-comparison/results/variant/Clair3-RNA/{bamtype}_{aligner}_{sample}.vcf.gz", 
    aligner=ONT_ALIGNER, sample=ONT_SAMPLE, bamtype=BAM)

deep_varaint = expand("results/variant/DeepVariant/{bamtype}_{aligner}_{sample}.vcf.gz", 
    aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, bamtype=BAM) + \
    expand("results/variant/DeepVariant/{bamtype}_{aligner}_{sample}.vcf.gz", 
    aligner=ONT_ALIGNER, sample=ONT_SAMPLE, bamtype=BAM)

longcallR = directory(expand("results/variant/longcallR/{bamtype}_{aligner}_{sample}", 
    aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, bamtype=BAM)) + \
    directory(expand("results/variant/longcallR/{bamtype}_{aligner}_{sample}", 
    aligner=ONT_ALIGNER, sample=ONT_SAMPLE, bamtype=BAM))

# vcfs =  expand("results/variant/{tool}/{bamtype}_{aligner}_{sample}.vcf.gz", 
#     aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, bamtype=BAM, tool=) +

# benchmarking
bed = expand("results/beds/{sample}/{aligner}_{coverage}.bed", aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, coverage=COVERAGE) + \
    expand("results/beds/{sample}/{aligner}_{coverage}.bed", aligner=ONT_ALIGNER, sample=ONT_SAMPLE, coverage=COVERAGE)

happy = expand("results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.summary.csv", tool=TOOL, aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, bamtype=BAM, coverage=COVERAGE) + \
    expand("results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.summary.csv", tool=TOOL,  aligner=ONT_ALIGNER, sample=ONT_SAMPLE, bamtype=BAM, coverage=COVERAGE)
# roc = expand("results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.roc.csv.gz",  tool=TOOL, aligner=ALIGNER, sample=SAMPLE, bamtype=BAM,  coverage=COVERAGE)
stratified = expand("results/stratified/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.extended.csv", tool=TOOL,  aligner=PACBIO_ALIGNER, sample=PACBIO_SAMPLE, bamtype=BAM, coverage=COVERAGE) + \
    expand("results/stratified/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.extended.csv", tool=TOOL,  aligner=ONT_ALIGNER, sample=ONT_SAMPLE, bamtype=BAM, coverage=COVERAGE)
res = {
    # "origin": origin,
    # "transformed": transformed,
    # "clair3-rna": clair3_rna,
    #"deep_variant": deep_varaint,
    #"longcallR": longcallR,
    "bed": bed,
    "happy": happy,
    # "roc": roc
    "stratified": stratified
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
rule bam2fq:
    priority: 99
    input:
        "data/bams/{sample}.flnc.bam"
    output:
        "results/align/genome/bam2fq/{sample}.fastq.gz"
    log:
        stderr="logs/align/genome/bam2fq/{sample}.stderr"
    threads: config["convert_threads"]
    params:
        threads=config["convert_threads"]
    shell:
        """
        samtools bam2fq -@ {params.threads} {input} | gzip > {output}
        """

def minimap2_preset(wildcards):
    if "IsoSeq" in wildcards.sample or "MasSeq" in wildcards.sample:
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
        "results/align/genome/origin/minimap2_{sample}.aligned.bam"
    params:
        align_map_bam_threads = config["align_map_bam_threads"],
        align_sort_bam_threads = config["align_sort_bam_threads"],
        align_sort_bam_memory_gb = config["align_sort_bam_memory_gb"],
        preset = minimap2_preset
    conda:
        "envs/align.yaml"
    log:
        stdout = "logs/minimap2/{sample}.out",
        stderr = "logs/minimap2/{sample}.err"
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
        reads = "results/align/genome/bam2fq/{sample}.fastq.gz",   
        genome = config["reference_genome"]
    output:
        bam = "results/align/genome/origin/pbmm2_{sample}.aligned.bam"
    params:
        preset = "CCS",  
        threads = config["align_map_bam_threads"]
    conda:
        "envs/align.yaml"
    log:
        stdout = "logs/pbmm2/{sample}.out",
        stderr = "logs/pbmm2/{sample}.err"
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



############################ Transformation ##################################
rule split_ncigar:
    priority: 97
    input:
        bam = "results/align/genome/origin/{aligner}_{sample}.aligned.bam",
        ref = config["reference_genome"],
        img = "/home/jiayiwang/tools/gatk-latest.simg"  
    output: 
        "results/align/genome/splitNC/{aligner}_{sample}.aligned.bam"
    log:
        "logs/split_ncigar_reads_{aligner}_{sample}.log"
    params:
        threads = 10
    shell:
        """
        singularity exec --bind /home/jiayiwang/miniconda3/envs/snakemake/bin/python:/usr/bin/python {input.img} /gatk/gatk --java-options "-Xmx16G -XX:+UseParallelGC -XX:ParallelGCThreads={params.threads} -Djava.util.concurrent.ForkJoinPool.common.parallelism={params.threads} -Djava.io.tmpdir=/home/jiayiwang/tmp1" SplitNCigarReads \
            -R {input.ref} \
            -I {input.bam} \
            -O {output} \
        > {log} 2>&1
        """ 

rule flag_correction:
    priority: 96
    input:
        bam = "results/align/genome/origin/{aligner}_{sample}.aligned.bam",
        sncr_bam = rules.split_ncigar.output
    output: 
        "results/align/genome/transformed/{aligner}_{sample}.aligned.bam"
    log:
        "logs/flag_correction_{aligner}_{sample}.log"
    params:
        threads = 5
    shell:
        """
        Rscript /home/jiayiwang/tools/lrRNAseqVariantCalling/tools/flagCorrection.r \
            {input.bam} \
            {input.sncr_bam} \
            {output} \
            {params.threads} \
        > {log} 2>&1
        samtools index {output}
        """

############################ Variant Calling ##################################

def get_clair3_platform(wildcards):
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
    elif "r9" in sample and "cdna" in sample:
        return "ont_r9_guppy_cdna"
    else:
        raise ValueError(f"Cannot determine platform for sample '{wildcards.sample}' and aligner '{wildcards.aligner}'")

rule run_clair3_rna:
    priority: 95
    input:
        bam="/home/jiayiwang/variant-comparison/results/align/genome/{bamtype}/{aligner}_{sample}.aligned.bam",
        ref=config["reference_genome"],
        img="/home/jiayiwang/tools/clair3-rna-latest.simg"
    output:
        dir=directory("/home/jiayiwang/variant-comparison/results/variant/Clair3-RNA/{bamtype}_{aligner}_{sample}"),
        file="/home/jiayiwang/variant-comparison/results/variant/Clair3-RNA/{bamtype}_{aligner}_{sample}.vcf.gz"
    log: 
        "logs/run_clair3_rna_{bamtype},{aligner},{sample}.log"
    params:
        threads = 10,
        model = get_clair3_platform
    shell:
        """
        if [ ! -f {input[0]}.bai ]; then
            samtools index {input[0]};
        fi
        singularity exec -B {input[0]},{input[1]} {input[2]} \
            /bin/bash -c "source /opt/conda/bin/activate /opt/conda/envs/clair3_rna && \
            /opt/bin/run_clair3_rna \
            --bam_fn {input[0]} \
            --ref_fn {input[1]} \
            --threads {params.threads} \
            --platform {params.model} \
            --tag_variant_using_readiportal \
            --output_dir {output[0]} \
            --conda_prefix /opt/conda/envs/clair3_rna" \
            > {log} 2>&1
        mv {output[0]}/output.vcf.gz {output[1]}
        """

# rule run_deep_variant:
#     priority: 95
#     input: 
#         bam="results/align/genome/{bamtype}/{aligner}_{sample}.aligned.bam",
#         ref=config["reference_genome"],
#         img="/home/jiayiwang/tools/deepvariant-1.8.0.simg"
#     output:
#         "results/variant/DeepVariant/{bamtype}_{aligner}_{sample}.vcf.gz",
#     log:
#         "logs/run_deep_variant_{bamtype},{aligner},{sample}.log"
#     params: 
#         threads = 8,
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
#             --num_shards {params.threads} \
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
        "results/align/genome/{bamtype}/{aligner}_{sample}.aligned.bam"
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

        if ! ~/tools/longcallR/target/release/longcallR \
            --bam-path {input} \
            --ref-path {params.ref} \
            --output {output.dir}/output \
            --preset {params.platform} \
            -t {params.threads}\
            --no-bam-output \
            --min-depth 0 \
            > {log} 2>&1; then
            echo "Error occurred, creating empty VCF file."
            echo "" | bgzip > {output.file}
            exit 0
        fi
        if [ -f {output.dir}/output.vcf ]; then
            bgzip {output.dir}/output.vcf
            mv {output.dir}/output.vcf.gz {output.file}
        else
            touch {output.file}
        fi
        """

############################ Benchmarking ##################################
rule generate_bed:
    priority: 93
    input: 
        bam="results/align/genome/transformed/{aligner}_{sample}.aligned.bam",
        gtf="/home/jiayiwang/pacbio/data/genome/gencode.v47.chr_patch_hapl_scaff.annotation.bed"
    output:
        "results/beds/{sample}/{aligner}_{coverage}.bed"
    params:
        threads = 8,
        img = "/home/jiayiwang/tools/mosdepth-0.3.3--h37c5b7d_2.simg",
        min_cov=lambda w: w.coverage.split('-')[0] if '-' in w.coverage else w.coverage[2:],
        max_cov=lambda w: w.coverage.split('-')[1] if '-' in w.coverage else None,
        confidence = lambda wildcards: f"data/truth/beds/{wildcards.sample.split('-')[0]}.bed",
    shell:
        """
        singularity exec {params.img} mosdepth \
            --threads {params.threads} \
            --fast-mode \
            results/beds/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage} \
            {input.bam}

        gzip -fdc results/beds/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.per-base.bed.gz \
        | awk -v min_cov={params.min_cov} -v max_cov={params.max_cov} \
            '$4 >= min_cov && (max_cov == "" || $4 <= max_cov)' \
        | bedtools merge -d 1 -c 4 -o mean -i - \
        > results/beds/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.coverage.bed

        bedtools intersect -a results/beds/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.coverage.bed \
            -b {input.gtf} > results/beds/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.tmp.bed

        if [ -f {params.confidence} ]; then
            bedtools intersect -a results/beds/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.tmp.bed \
                -b {params.confidence} > {output}
        else
            mv results/beds/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.tmp.bed {output}
        fi
        """



rule happy_benchmark:
    priority: 92
    input:
        vcf="results/variant/{tool}/{bamtype}_{aligner}_{sample}.vcf.gz",
        bed="results/beds/{sample}/{aligner}_{coverage}.bed"
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
        export HGREF=/home/jiayiwang/pacbio/data/genome/GRCh38.p14.genome.fa
        
        # Try to run the full pipeline
        if (
           
            if ! bcftools index -t -f "{input.vcf}"; then
                echo "VCF is not sorted. Sorting and reindexing..."
                cp {input.vcf} {input.vcf}.unsorted
                bcftools sort -Oz -o {input.vcf} {input.vcf}.unsorted
                bcftools index -t -f {input.vcf}
                rm {input.vcf}.unsorted
            fi
            singularity exec \
                --bind /home/jiayiwang/variant-comparison:/home/jiayiwang/variant-comparison \
                --bind /home/jiayiwang/pacbio:/home/jiayiwang/pacbio \
                /home/jiayiwang/tools/hap.py.simg /opt/hap.py/bin/hap.py \
                "{params.truth}" \
                "{input.vcf}" \
                -r "{params.ref}" \
                -o "{params.outdir}" \
                --engine=vcfeval \
                --target-regions="{input.bed}" \
                --threads="{params.threads}"
            
            # Process the output if hap.py succeeded
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
        )
        then
            echo "Pipeline completed successfully"
        else
            echo "Pipeline failed - keeping empty output file"
        fi
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