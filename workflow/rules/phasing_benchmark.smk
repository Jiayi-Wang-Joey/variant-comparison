TYPE = ["all", "SNP"]
N_VARIANTS = ["TOTAL"]

rule subset_variant:
    input:
        vcf = "results/phase/{phaser}/{tool}_origin_minimap2_{sample}.vcf.gz"
    output:
        vcf = "results/phase_typed/{phaser}/{tool}_origin_minimap2_{sample}_{type}.vcf.gz"
    log:
        "logs/subset_variant_{phaser}_{tool}_{sample}_{type}.log"
    conda:
        "../envs/phasing_benchmark.yaml"
    shell:
        """

        if [ "{wildcards.type}" = "SNP" ]; then
            bcftools view -v snps -O z -o "{output.vcf}" "{input.vcf}" >> "{log}" 2>&1
        else
            cp "{input.vcf}" "{output.vcf}" >> "{log}" 2>&1
        fi

        tabix -f -p vcf "{output.vcf}" >> "{log}" 2>&1
        """
    
rule longcallR_variant:
    input:
        vcf = "results/variant/longcallR/origin_minimap2_{sample}.vcf.gz"
    output:
        vcf = "results/phase_typed/longcallR/longcallR_origin_minimap2_{sample}_{type}.vcf.gz"
    log:
        "logs/subset_variant_longcallR_longcallR_{sample}_{type}.log"
    conda:
        "../envs/phasing_benchmark.yaml"
    shell:
        """
        mkdir -p tmp

        SORTED=tmp/longcallR_origin_minimap2_{wildcards.sample}.sorted.vcf.gz
        bcftools sort -Oz -o "$SORTED" "{input.vcf}" >> "{log}" 2>&1

        if [ "{wildcards.type}" = "SNP" ]; then
            bcftools view -v snps -O z -o "{output.vcf}" "$SORTED" >> "{log}" 2>&1
        else
            cp "$SORTED" "{output.vcf}" >> "{log}" 2>&1
        fi

        tabix -f -p vcf "{output.vcf}" >> "{log}" 2>&1
        rm -f "$SORTED"
        """


rule isolaser_variant:
    input:
        vcf = "results/variant/isoLASER_filtered/origin_minimap2_{sample}.vcf.gz"
    output:
        vcf = "results/phase_typed/isoLASER/isoLASER_origin_minimap2_{sample}_{type}.vcf.gz"
    log:
        "logs/subset_variant_isoLASER_isoLASER_{sample}_{type}.log"
    conda:
        "../envs/phasing_benchmark.yaml"
    params:
        script = "workflow/scripts/convert_isolaser_vcf.py",
        tmp = "tmp/isoLASER_origin_minimap2_{sample}_converted.vcf.gz"
    shell:
        """
        mkdir -p tmp

        # Convert isoLASER gVCF: remove NON_REF, deduplicate positions, add PS tag
        python3 {params.script} {input.vcf} {wildcards.sample} {params.tmp} >> {log} 2>&1
        tabix -f -p vcf {params.tmp} >> {log} 2>&1

        bcftools view -v snps -Oz -o {output.vcf} {params.tmp} >> {log} 2>&1

        tabix -f -p vcf {output.vcf} >> {log} 2>&1
        rm -f {params.tmp} {params.tmp}.tbi
        """


def get_bed_phasing(wildcards):
    sample = wildcards.sample
    tool = wildcards.tool
    if "IsoSeq" in sample and tool == "DeepVariant":
        status = "preprocessed"
    else:
        status = "raw"
    return f"results/beds/{status}/{sample}_minimap2_5.bed"


rule compute_switch_error:
    input:
        truth = lambda wc: f"data/truth/{wc.sample[:5]}/phase.vcf.gz",
        query = "results/phase_typed/{phaser}/{tool}_origin_minimap2_{sample}_{type}.vcf.gz",
        bed = get_bed_phasing
    output:
        tsv = "results/switch_error/raw/{phaser}_{tool}_origin_minimap2_{sample}_{type}.tsv"
    log:
        "logs/switch_error_{phaser}_{tool}_{sample}_{type}.log"
    conda:
        "../envs/whatshap.yaml"
    shell:
        """
        tabix -f -p vcf "{input.truth}" >> "{log}" 2>&1
        tabix -f -p vcf "{input.query}" >> "{log}" 2>&1

        TRUTH_R=tmp/{wildcards.phaser}_{wildcards.tool}_{wildcards.sample}_{wildcards.type}.truth.restricted.vcf.gz
        QUERY_R=tmp/{wildcards.phaser}_{wildcards.tool}_{wildcards.sample}_{wildcards.type}.query.restricted.vcf.gz

        bcftools view -R "{input.bed}" -Oz -o "$TRUTH_R" "{input.truth}" >> "{log}" 2>&1
        bcftools view -R "{input.bed}" -Oz -o "$QUERY_R" "{input.query}" >> "{log}" 2>&1
        tabix -f -p vcf "$TRUTH_R" >> "{log}" 2>&1
        tabix -f -p vcf "$QUERY_R" >> "{log}" 2>&1

        whatshap compare --ignore-sample-name \
            --tsv-pairwise "{output.tsv}" \
            "$TRUTH_R" "$QUERY_R" >> "{log}" 2>&1

        rm -f "$TRUTH_R" "$TRUTH_R.tbi" "$QUERY_R" "$QUERY_R.tbi"
        """


rule summarize_switch_error:
    input:
        tsv = "results/switch_error/raw/{phaser}_{tool}_origin_minimap2_{sample}_{type}.tsv",
        vcf = "results/phase_typed/{phaser}/{tool}_origin_minimap2_{sample}_{type}.vcf.gz"
    output:
        tsv = "results/switch_error/summary/{phaser}_{tool}_origin_minimap2_{sample}_{type}.tsv"
    log:
        "logs/summarize_switch_error_{phaser}_{tool}_{sample}_{type}.log"
    conda:
        "../envs/phasing_benchmark.yaml"
    shell:
        """
        set -euo pipefail

        n_actual=$(bcftools view -H "{input.vcf}" 2> "{log}" | wc -l)
        n_phased=$(bcftools view -H -i 'GT~"|"' "{input.vcf}" | wc -l)
        echo "Sample: {wildcards.sample}, Caller: {wildcards.tool}, Phaser: {wildcards.phaser}, Type: {wildcards.type}" > "{log}"
        awk -v sample="{wildcards.sample}" \
            -v caller="{wildcards.tool}" \
            -v phaser="{wildcards.phaser}" \
            -v type="{wildcards.type}" \
            -v n_actual="$n_actual" \
            -v n_phased="$n_phased" \
            'BEGIN {{FS=OFS="\\t"}}
            NR==1 {{
                print "sample","caller","phaser","type","n_actual","n_phased",$0
            }}
            NR>1 {{
                print sample,caller,phaser,type,n_actual,n_phased,$0
            }}' \
            {input.tsv} > {output.tsv}
        """


rule phase_block_stats:
    input:
        vcf = "results/phase_typed/{phaser}/{tool}_origin_minimap2_{sample}_{type}.vcf.gz",
        bed = get_bed_phasing
    output:
        stats  = "results/phase_blocks/raw/{phaser}_{tool}_origin_minimap2_{sample}_{type}.stats.tsv",
        blocks = "results/phase_blocks/raw/{phaser}_{tool}_origin_minimap2_{sample}_{type}.blocks.tsv",
        gtf    = "results/phase_blocks/raw/{phaser}_{tool}_origin_minimap2_{sample}_{type}.blocks.gtf"
    log:
        "logs/phase_blocks_{phaser}_{tool}_{sample}_{type}.log"
    conda:
        "../envs/whatshap.yaml"
    shell:
        """
        tabix -f -p vcf "{input.vcf}" >> "{log}" 2>&1

        VCF_R=tmp/{wildcards.phaser}_{wildcards.tool}_{wildcards.sample}_{wildcards.type}.phaseblocks.restricted.vcf.gz

        bcftools view -R "{input.bed}" -Oz -o "$VCF_R" "{input.vcf}" >> "{log}" 2>&1
        tabix -f -p vcf "$VCF_R" >> "{log}" 2>&1

        whatshap stats \
            --tsv "{output.stats}" \
            --block-list "{output.blocks}" \
            --gtf "{output.gtf}" \
            "$VCF_R" >> "{log}" 2>&1

        rm -f "$VCF_R" "$VCF_R.tbi"
        """


rule summarize_phase_blocks:
    input:
        stats  = rules.phase_block_stats.output.stats,
        blocks = rules.phase_block_stats.output.blocks,
        gtf    = rules.phase_block_stats.output.gtf
    output:
        tsv = "results/phase_blocks/summary/{phaser}_{tool}_origin_minimap2_{sample}_{type}.tsv"
    log:
        "logs/summarize_phase_blocks_{phaser}_{tool}_{sample}_{type}.log"
    conda:
        "../envs/phasing_benchmark.yaml"
    params:
        script = "workflow/scripts/summarize_phase_blocks.py"
    shell:
        """
        python3 {params.script} \
            --blocks "{input.blocks}" \
            --gtf "{input.gtf}" \
            --stats "{input.stats}" \
            --sample "{wildcards.sample}" \
            --caller "{wildcards.tool}" \
            --phaser "{wildcards.phaser}" \
            --type "{wildcards.type}" \
            --output "{output.tsv}" > "{log}" 2>&1
        """

