rule generate_bed:
    priority: 93
    input: 
        bam="results/align/{status}/origin/{aligner}_{sample}.aligned.bam",
        gtf=config["gtf"]
    output:
        "results/beds/{status}/{sample}_{aligner}_{coverage}.bed"
    params:
        threads = 8,
        img = "/home/jiayiwang/tools/mosdepth-0.3.3--h37c5b7d_2.simg",
        min_cov=lambda w: w.coverage,
        confidence = lambda wildcards: f"data/truth/beds/{wildcards.sample.split('-')[0].rstrip('abc')}.bed",
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
        confidence = lambda wildcards: f"data/truth/beds/{wildcards.sample.split('-')[0].rstrip('abc')}.bed"
    output:
        smy="results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.summary.csv",
        roc="results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.roc.csv.gz",
        vcf="results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.vcf.gz"
    params:
        truth = lambda wildcards: f"data/truth/{wildcards.sample.split('-')[0].rstrip('abc')}/truth.vcf.gz",
        ref = config["reference_genome"],
        outdir = "results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}",
        threads = 10 
    conda:
        "../envs/happy.yaml"
    log:
        "logs/happy_{tool}_{bamtype}_{aligner}_{sample}_{coverage}.log"
    shell:
        """
            export HGREF={params.ref}
        
            if ! bcftools index -t -f "{input.vcf}" >/dev/null 2>&1; then
                echo "VCF is not sorted or not indexed properly. Sorting and indexing..."
                cp "{input.vcf}" "{input.vcf}.unsorted"
                bcftools sort -Oz -o "{input.vcf}.tmp.gz" "{input.vcf}.unsorted"
                mv "{input.vcf}.tmp.gz" "{input.vcf}"
                bcftools index -t -f "{input.vcf}"
                rm "{input.vcf}.unsorted"
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
        stratification="/home/jiayiwang/variant-comparison/data/regions/GRCh38-all-stratifications.tsv"
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

rule annotate_vcfs:
    input:
        vcf = "results/happy/{tool}_{bamtype}_{aligner}_{sample}_5.vcf.gz"
    output:
        "results/SnpEff/{tool}_{bamtype}_{aligner}_{sample}_5.vcf.gz"
    log: 
        "logs/SnpEff_{tool}_{bamtype}_{aligner}_{sample}_5.vcf.gz"
    shell:
        """
        java -Xmx8g -jar ~/snpEff/snpEff.jar GRCh38.86 {input.vcf} | bgzip -c > {output}
        """