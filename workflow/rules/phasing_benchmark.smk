TYPE = ["all", "SNP"]
N_VARIANTS = ["TOTAL"]

rule subset_variant:
    input:
        vcf = "results/phase/{phaser}/{tool}_origin_minimap2_{sample}.vcf.gz"
    output:
        vcf = "results/phase_typed/{phaser}/{tool}_origin_minimap2_{sample}_{type}.vcf.gz"
    log:
        "logs/subset_variant_{phaser}_{tool}_{sample}_{type}.log"
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
    shell:
        """
        if [ "{wildcards.type}" = "SNP" ]; then
            bcftools view -v snps -O z -o "{output.vcf}" "{input.vcf}" >> "{log}" 2>&1
        else
            cp "{input.vcf}" "{output.vcf}" >> "{log}" 2>&1
        fi

        tabix -f -p vcf "{output.vcf}" >> "{log}" 2>&1
        """


rule compute_switch_error:
    input:
        truth = lambda wc: f"data/truth/{wc.sample[:5]}/phase.vcf.gz",
        query = "results/phase_typed/{phaser}/{tool}_origin_minimap2_{sample}_{type}.vcf.gz"
    output:
        tsv = "results/switch_error/raw/{phaser}_{tool}_origin_minimap2_{sample}_{type}.tsv"
    log:
        "logs/switch_error_{phaser}_{tool}_{sample}_{type}.log"
    shell:
        """

        tabix -f -p vcf "{input.truth}" >> "{log}" 2>&1
        tabix -f -p vcf "{input.query}" >> "{log}" 2>&1

        whatshap compare --ignore-sample-name \
            --tsv-pairwise "{output.tsv}" \
            "{input.truth}" "{input.query}" >> "{log}" 2>&1

        """


rule summarize_switch_error:
    input:
        tsv = "results/switch_error/raw/{phaser}_{tool}_origin_minimap2_{sample}_{type}.tsv",
        vcf = "results/phase_typed/{phaser}/{tool}_origin_minimap2_{sample}_{type}.vcf.gz"
    output:
        tsv = "results/switch_error/summary/{phaser}_{tool}_origin_minimap2_{sample}_{type}.tsv"
    log:
        "logs/summarize_switch_error_{phaser}_{tool}_{sample}_{type}.log"
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

