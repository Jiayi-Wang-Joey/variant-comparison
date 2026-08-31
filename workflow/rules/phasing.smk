def get_bam(wildcards):
    sample = wildcards.sample
    tool = wildcards.tool
    aligner = "minimap2"

    if "IsoSeq" in sample and tool == "DeepVariant":
        status = "preprocessed"
    else:
        status = "raw"

    return f"results/align/{status}/origin/{aligner}_{sample}.aligned.bam"


# rule filter_vcf:
#     input:
#         vcf = "results/variant/{tool}/origin_minimap2_{sample}.vcf.gz"
#     output:
#         vcf = "results/variant/{tool}/filtered_{sample}.vcf.gz"
#     shell:
#         "bcftools view -f PASS {input.vcf} -Oz -o {output.vcf} && bcftools index {output.vcf}"

rule filter_vcf:
    input:
        vcf = "results/variant/{tool}/origin_minimap2_{sample}.vcf.gz",
        ref = config["reference_genome"]
    output:
        vcf = "results/variant/{tool}/filtered_{sample}.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools view -f PASS {input.vcf} \
        | bcftools norm -f {input.ref} -c s -Oz -o {output.vcf}
        
        bcftools index -f {output.vcf}
        """

rule WhatsHap:
    input:
        vcf = rules.filter_vcf.output,
        bam = get_bam,
        ref = config["reference_genome"]
    output:
        "results/phase/WhatsHap/{tool}_origin_minimap2_{sample}.vcf.gz"
    log:
        "logs/phase_WhatsHap_{tool}_origin_minimap2_{sample}.log"
    conda:
        "../envs/whatshap.yaml"
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
        vcf = rules.filter_vcf.output,
        bam = get_bam,
        ref = config["reference_genome"]
    output:
        "results/phase/longphase/{tool}_origin_minimap2_{sample}.vcf.gz"
    threads: 8
    params:
        platform = lambda wildcards: get_platform_longphase(wildcards)
    log:
        "logs/phase_longphase_{tool}_origin_minimap2_{sample}.log"
    conda:
        "../envs/bcftools.yaml"
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

# rule longcallR_phase:
#     input:
#         vcf = rules.filter_vcf.output,
#         bam = get_bam,
#         ref = config["reference_genome"]
#     output:
#         "results/phase/longcallR_phase/{tool}_origin_minimap2_{sample}.vcf.gz"
#     params:
#         path = directory("results/phase/longcallR_phase/{tool}_origin_minimap2_{sample}"),
#         platform = get_longcallR_platform,
#         threads = 8
#     log: 
#         "logs/phase_longcallR_{tool}_origin_minimap2_{sample}.log"
#     shell:
#         """
#         mkdir -p {params.path}
#         ~/tools/longcallR/target/release/longcallR --bam-path {input.bam} \
#         --input-vcf {input.vcf} \
#         --ref-path {input.ref} \
#         --preset {params.platform} \
#         --no-bam-output \
#         -t {params.threads} \
#         --output {params.path}/output \
#         > {log} 2>&1

#         bgzip {params.path}/output.vcf
#         mv {params.path}/output.vcf.gz {output}
#         """


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
        vcf = rules.filter_vcf.output,
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
        set -euo pipefail

        mkdir -p tmp/{wildcards.tool}/

        VCF_BASE=tmp/{wildcards.tool}/origin_minimap2_{wildcards.sample}.hapcut2.vcf

        # 1. Normalize and clean VCF for HapCUT2.
        #    extractHAIRS asserts the FORMAT column starts with "GT". isoLASER places other
        #    fields before GT, so for isoLASER we strip FORMAT down to GT-only (plain text).
        #    HAPCUT2 must be given the EXACT SAME file as extractHAIRS, so a single file is
        #    used for both tools in each branch.
        if [ "{wildcards.tool}" = "isoLASER" ]; then
            bcftools view -g ^miss {input.vcf} \
              | bcftools norm -m-any \
              | bcftools view -Ov \
                  -i 'GT="0/0" || GT="0/1" || GT="1/1" || GT="1/2" || GT="2/2" || GT="0|1" || GT="1|0"' \
              | awk 'BEGIN{{OFS="\t"}}
                    /^##FORMAT=<ID=GT,/{{print; next}}
                    /^##FORMAT/{{next}}
                    /^#/{{print; next}}
                    {{
                      n=split($9,F,":");
                      gtidx=0;
                      for(i=1;i<=n;i++) if(F[i]=="GT"){{gtidx=i; break}}
                      if(gtidx==0) next;
                      split($10,S,":");
                      $9="GT"; $10=S[gtidx]; print
                    }}' \
              > "$VCF_BASE"
        else
            # bcftools annotate -x errors out (rather than no-op) if none of the
            # requested tags exist in the header, so only request tags that are present.
            TAGS=""
            if bcftools view -h {input.vcf} | grep -q "^##FORMAT=<ID=VAF,"; then TAGS="FORMAT/VAF"; fi
            if bcftools view -h {input.vcf} | grep -q "^##FORMAT=<ID=PL,"; then
                if [ -n "$TAGS" ]; then TAGS="$TAGS,FORMAT/PL"; else TAGS="FORMAT/PL"; fi
            fi

            if [ -n "$TAGS" ]; then
                bcftools view -g ^miss {input.vcf} \
                  | bcftools norm -m-any \
                  | bcftools annotate -x "$TAGS" \
                  | bcftools view -i 'GT="0/0" || GT="0/1" || GT="1/1" || GT="1/2" || GT="2/2"' \
                  -Ov -o "$VCF_BASE"
            else
                bcftools view -g ^miss {input.vcf} \
                  | bcftools norm -m-any \
                  | bcftools view -i 'GT="0/0" || GT="0/1" || GT="1/1" || GT="1/2" || GT="2/2"' \
                  -Ov -o "$VCF_BASE"
            fi
        fi

        # 2. Extract fragments
        extractHAIRS {params.platform} 1 \
            --bam {input.bam} \
            --VCF "$VCF_BASE" \
            --out tmp/{wildcards.sample}.fragment_file \
            --ref {input.ref} >> {log} 2>&1

        # 3. Run HapCUT2 phasing. --outvcf 1 is REQUIRED for this build to emit
        #    <output>.phased.vcf (default is 0 = no VCF output). Must use the same --VCF
        #    file that was given to extractHAIRS.
        HAPCUT2 --fragments tmp/{wildcards.sample}.fragment_file \
            --VCF "$VCF_BASE" \
            --outvcf 1 \
            --output results/phase/HapCUT2/{wildcards.tool}_origin_minimap2_{wildcards.sample} >> {log} 2>&1

        # 4. Index phased VCF (HapCUT2 writes .phased.vcf; glob covers either case)
        PHASED_VCF=$(ls results/phase/HapCUT2/{wildcards.tool}_origin_minimap2_{wildcards.sample}.phased.[Vv][Cc][Ff] 2>/dev/null | head -1)
        bgzip -f "$PHASED_VCF"
        mv "${{PHASED_VCF}}.gz" {output.vcf}
        tabix -p vcf {output.vcf} >> {log} 2>&1

        # 5. Cleanup
        rm -f "$VCF_BASE" \
              tmp/{wildcards.sample}.fragment_file
        """

rule HiPhase:
    input:
        vcf = rules.filter_vcf.output,
        bam = get_bam,
        ref = config["reference_genome"]
    output:
        vcf = "results/phase/HiPhase/{tool}_origin_minimap2_{sample}.vcf.gz"
    log:
        "logs/phase_HiPhase_{tool}_origin_minimap2_{sample}.log"
    threads: 8
    conda:
        "../envs/phase.yaml"
    shell:
        """
        hiphase \
            --bam {input.bam} \
            --vcf {input.vcf} \
            --preset rna \
            --output-vcf {output.vcf} \
            --reference {input.ref} \
            --ignore-read-groups \
            --phase-singletons \
            --threads {threads} &>> {log}
        """

# def get_margin_params(wildcards):
#     sample = wildcards.sample
#     if sample in PACBIO_SAMPLE:
#         return "/opt/margin_dir/params/phase/allParams.phase_vcf.pb-hifi.json"
#     elif sample in ONT_SAMPLE:
#         return "/opt/margin_dir/params/phase/allParams.phase_vcf.ont.json"
#     else:
#         raise ValueError(f"Unknown platform for sample: {sample}. Check PACBIO_SAMPLE/ONT_SAMPLE lists.")

# rule margin:
#     input:
#         vcf = rules.filter_vcf.output.vcf,
#         bam = get_bam,
#         ref = config["reference_genome"]
#     output:
#         vcf = "results/phase/Margin/{tool}_origin_minimap2_{sample}.vcf.gz"
#     log:
#         "logs/phase_Margin_{tool}_origin_minimap2_{sample}.log"
#     threads: 8
#     container:
#         "docker://kishwars/pepper_deepvariant:r0.8"
#     params:
#         path = "results/phase/Margin/{tool}_origin_minimap2_{sample}",
#         vcf = "results/phase/Margin/{tool}_origin_minimap2_{sample}.phased.vcf",
#         margin_json = get_margin_params
#     shell:
#         """
#         margin phase \
#             {input.bam} \
#             {input.ref} \
#             {input.vcf} \
#             {params.margin_json} \
#             -t {threads} \
#             --skipHaplotypeBAM \
#             -o {params.path} > {log} 2>&1

#         bgzip -f {params.vcf} &>> {log}
#         mv {params.vcf}.gz {output.vcf}
#         """