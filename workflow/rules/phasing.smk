def get_bam(wildcards):
    sample = wildcards.sample
    tool = wildcards.tool
    aligner = "minimap2"

    if "IsoSeq" in sample and tool == "DeepVariant":
        status = "preprocessed"
    else:
        status = "raw"

    return f"results/align/{status}/origin/{aligner}_{sample}.aligned.bam"

rule WhatsHap:
    input:
        vcf = "results/variant/{tool}/origin_minimap2_{sample}.vcf.gz",
        bam = get_bam,
        ref = config["reference_genome"]
    output:
        "results/phase/WhatsHap/{tool}_origin_minimap2_{sample}.vcf.gz"
    log:
        "logs/phase_WhatsHap_{tool}_origin_minimap2_{sample}.log"
    shell:
        """
        whatshap phase --ignore-read-groups \
            -o {output} \
            --reference={input.ref} \
            {input.vcf} {input.bam} > {log} 2>&1
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
        vcf = "results/variant/{tool}/origin_minimap2_{sample}.vcf.gz",
        bam = get_bam,
        ref = config["reference_genome"]
    output:
        "results/phase/longphase/{tool}_origin_minimap2_{sample}.vcf.gz"
    threads: 8
    params:
        platform = lambda wildcards: get_platform_longphase(wildcards)
    log:
        "logs/phase_longphase_{tool}_origin_minimap2_{sample}.log"
    shell:
        """
        /home/jiayiwang/tools/longphase_linux-x64 phase \
            -s {input.vcf} \
            -b {input.bam} \
            -r {input.ref} \
            -t {threads} \
            -o results/phase/longphase/{wildcards.tool}_origin_minimap2_{wildcards.sample} \
            --indels \
            {params.platform} \
            > {log} 2>&1

        bgzip results/phase/longphase/{wildcards.tool}_origin_minimap2_{wildcards.sample}.vcf
        tabix -p vcf results/phase/longphase/{wildcards.tool}_origin_minimap2_{wildcards.sample}.vcf.gz
        """

rule longcallR_phase:
    input:
        vcf = "results/variant/{tool}/origin_minimap2_{sample}.vcf.gz",
        bam = get_bam,
        ref = config["reference_genome"]
    output:
        "results/phase/longcallR_phase/{tool}_origin_minimap2_{sample}.vcf.gz"
    params:
        path = directory("results/phase/longcallR_phase/{tool}_origin_minimap2_{sample}"),
        platform = get_longcallR_platform,
        threads = 8
    log: 
        "logs/phase_longcallR_{tool}_origin_minimap2_{sample}.log"
    shell:
        """
        mkdir -p {params.path}
        ~/tools/longcallR/target/release/longcallR --bam-path {input.bam} \
        --input-vcf {input.vcf} \
        --ref-path {input.ref} \
        --preset {params.platform} \
        --no-bam-output \
        -t {params.threads} \
        --output {params.path}/output \
        > {log} 2>&1

        bgzip {params.path}/output.vcf
        mv {params.path}/output.vcf.gz {output}
        """


def get_platform_hapcut2(wildcards):
    sample = wildcards.sample
    if sample in PACBIO_SAMPLE:
        return "--pacbio"
    elif sample in ONT_SAMPLE:
        return "--ont"
    else:
        raise ValueError(f"Unknown platform for sample: {sample}")


rule hapcut2:
    input:
        vcf = "results/variant/{tool}/origin_minimap2_{sample}.vcf.gz",
        bam = get_bam,
        ref = config["reference_genome"]
    output:
        vcf = "results/phase/HapCUT2/{tool}_origin_minimap2_{sample}.vcf.gz"
    log:
        "logs/phase_hapcut2_{tool}_origin_minimap2_{sample}.log"
    threads: 8
    conda:
        "../envs/phase.yaml"
    params:
        platform = lambda wildcards: get_platform_hapcut2(wildcards)
    shell:
        """
        mkdir -p tmp/{wildcards.tool}/

        # 1. Normalize and clean VCF for HapCUT2
        bcftools view -g ^miss {input.vcf} \
          | bcftools norm -m-any \
          | bcftools annotate -x FORMAT/VAF,FORMAT/PL \
          | bcftools view -i 'GT="0/0" || GT="0/1" || GT="1/1" || GT="1/2" || GT="2/2"' \
          -Oz -o tmp/{wildcards.tool}/origin_minimap2_{wildcards.sample}.hapcut2.vcf

        # 2. Extract fragments
        extractHAIRS {params.platform} 1 \
            --bam {input.bam} \
            --VCF tmp/{wildcards.tool}/origin_minimap2_{wildcards.sample}.hapcut2.vcf \
            --out tmp/{wildcards.sample}.fragment_file \
            --ref {input.ref} >> {log} 2>&1

        # 3. Run HapCUT2 phasing
        HAPCUT2 --fragments tmp/{wildcards.sample}.fragment_file \
            --VCF tmp/{wildcards.tool}/origin_minimap2_{wildcards.sample}.hapcut2.vcf \
            --output results/phase/HapCUT2/{wildcards.tool}_origin_minimap2_{wildcards.sample} >> {log} 2>&1

        # 4. Index phased VCF
        bgzip -f results/phase/HapCUT2/{wildcards.tool}_origin_minimap2_{wildcards.sample}.phased.VCF 
        mv results/phase/HapCUT2/{wildcards.tool}_origin_minimap2_{wildcards.sample}.phased.VCF.gz {output.vcf}
        tabix -p vcf {output.vcf} >> {log} 2>&1

        # 5. Cleanup
        rm -f tmp/{wildcards.tool}/origin_minimap2_{wildcards.sample}.hapcut2.vcf.gz* \
              tmp/{wildcards.sample}.fragment_file
        """

rule switch_error:
    input:
        vcf = "results/phase/{phaser}/{tool}_origin_minimap2_{sample}.vcf.gz"
    output:
        bed = "results/switch_error/{phaser}_{tool}_origin_minimap2_{sample}.bed",
        tsv = "results/switch_error/{phaser}_{tool}_origin_minimap2_{sample}.tsv"
    log:
        "logs/phase_switch_error_{phaser}_{tool}_origin_minimap2_{sample}.log"
    params:
        truth = lambda wildcards: f"data/truth/{wildcards.sample[:5]}/phase.vcf.gz"
    shell:
        """
        # Index truth VCF if missing
        if [ ! -f {params.truth}.tbi ]; then
            echo "Index file for truth VCF missing, creating with tabix..."
            tabix -p vcf {params.truth}
        fi

        # Index input VCF if missing
        if [ ! -f {input.vcf}.tbi ]; then
            echo "Index file for input VCF missing, sorting and indexing..."
            bcftools sort -O z -o {input.vcf}.sorted.vcf.gz {input.vcf}
            mv {input.vcf}.sorted.vcf.gz {input.vcf}
            tabix -p vcf {input.vcf}
        fi

        whatshap compare --ignore-sample-name \
            --switch-error-bed {output.bed} \
            --tsv-pairwise {output.tsv}.tmp \
            {params.truth} {input.vcf} > {log} 2>&1

        awk -v sample="{wildcards.sample}" \
            -v bamtype="origin" \
            -v caller="{wildcards.tool}" \
            -v phaser="{wildcards.phaser}" \
            -v aligner="minimap2" \
            'BEGIN {{FS=OFS="\\t"}} NR==1 {{print "sample","bamtype","caller","phaser","aligner",$0}} NR>1 {{print sample,bamtype,caller,phaser,aligner,$0}}' \
            {output.tsv}.tmp > {output.tsv}

        rm {output.tsv}.tmp
        """
