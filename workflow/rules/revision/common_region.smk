# hap.py on the common region (standard BED ∩ isoLASER effective-depth BED)


rule isolaser_anno_bed:
    priority: 50
    input:
        bam="results/variant/isoLASER/realign/origin/{aligner}_{sample}.anno.bam",
        gtf=config["gtf"]
    output:
        "results/beds/isolaser_anno/{sample}_{aligner}_{coverage}.bed"
    conda:
        "../../envs/bcftools.yaml"
    params:
        threads=8,
        img="/home/jiayiwang/tools/mosdepth-0.3.3--h37c5b7d_2.simg",
        prefix="results/beds/isolaser_anno/{sample}/{aligner}_{coverage}",
        confidence=lambda wc: f"data/truth/beds/{wc.sample.split('-')[0].rstrip('abc')}.bed",
    shell:
        """
        if [ ! -f {input.bam}.bai ] || [ {input.bam} -nt {input.bam}.bai ]; then
            samtools index {input.bam}
        fi
        mkdir -p $(dirname {params.prefix})
        singularity exec \
            --bind $(pwd -P):$(pwd -P) \
            --pwd $(pwd -P) \
            {params.img} mosdepth \
            --threads {params.threads} \
            --mapq 20 \
            {params.prefix} \
            {input.bam}

        gzip -fdc {params.prefix}.per-base.bed.gz \
            | awk -v min_cov={wildcards.coverage} '$4 >= min_cov' \
            | bedtools merge -d 1 -c 4 -o mean -i - \
            | bedtools intersect -a - -b <(grep -v "^#" {input.gtf} | awk '$3 == "exon"') \
            | bedtools intersect -a - -b {params.confidence} \
            > {output}
        rm -f {params.prefix}.per-base.bed.gz {params.prefix}.per-base.bed.gz.csi
        """


rule common_bed:
    priority: 50
    input:
        "results/beds/raw/{sample}_{aligner}_{coverage}.bed",
        "results/beds/isolaser_anno/{sample}_{aligner}_{coverage}.bed"
    output:
        "results/beds/common/{sample}_{aligner}_{coverage}.bed"
    conda:
        "../../envs/bcftools.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        # merge first: input BEDs repeat intervals per exon record
        sort -k1,1 -k2,2n {input[0]} | bedtools merge -i - > {output}.a
        sort -k1,1 -k2,2n {input[1]} | bedtools merge -i - > {output}.b
        bedtools intersect -a {output}.a -b {output}.b | sort -k1,1 -k2,2n | bedtools merge -i - > {output}
        rm -f {output}.a {output}.b
        """


rule happy_common:
    priority: 0
    retries: 2
    input:
        vcf=get_variant_vcf,
        bed="results/beds/common/{sample}_{aligner}_{coverage}.bed",
        confidence=lambda wc: f"data/truth/beds/{wc.sample.split('-')[0].rstrip('abc')}.bed"
    output:
        smy="results/happy_common/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.summary.csv",
        roc="results/happy_common/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.roc.csv.gz",
        vcf="results/happy_common/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.vcf.gz"
    wildcard_constraints:
        tool="[^/_]+"
    conda:
        "../../envs/bcftools.yaml"
    params:
        truth=lambda wc: f"data/truth/{wc.sample.split('-')[0].rstrip('abc')}/truth.vcf.gz",
        ref=config["reference_genome"],
        outdir="results/happy_common/{tool}_{bamtype}_{aligner}_{sample}_{coverage}",
        threads=10
    log:
        "logs/happy_common_{tool}_{bamtype}_{aligner}_{sample}_{coverage}.log"
    shell:
        """
        export HGREF={params.ref}
        mkdir -p results/happy_common tmp

        JOBTAG=common_{wildcards.tool}_{wildcards.bamtype}_{wildcards.aligner}_{wildcards.sample}_{wildcards.coverage}
        TMP=tmp/${{JOBTAG}}.vcf.gz

        # GATK QUAL>=30, others PASS
        if [ "{wildcards.tool}" = "GATK" ]; then
            bcftools view -i 'QUAL>=30' {input.vcf} | bcftools sort -Oz -o $TMP
        else
            bcftools view -f PASS {input.vcf} | bcftools sort -Oz -o $TMP
        fi
        bcftools index -t -f $TMP

        export TMPDIR=/data/jiayiwang/tmp
        mkdir -p $TMPDIR
        singularity exec \
            --bind $(pwd -P):$(pwd -P) \
            --bind $TMPDIR:$TMPDIR \
            --pwd $(pwd -P) \
            /home/jiayiwang/tools/hap.py.simg /opt/hap.py/bin/hap.py \
            {params.truth} \
            $TMP \
            -r {params.ref} \
            -f {input.confidence} \
            -T {input.bed} \
            -o {params.outdir} \
            --engine=vcfeval \
            --scratch-prefix $TMPDIR \
            --threads={params.threads} > {log} 2>&1

        awk -v tool="{wildcards.tool}" -v bamtype="{wildcards.bamtype}" -v aligner="{wildcards.aligner}" \
            -v sample="{wildcards.sample}" -v coverage="{wildcards.coverage}" \
        'BEGIN {{FS=OFS=","}}
            NR==1 {{print $0, "tool", "bamtype", "aligner", "sample", "coverage"}}
            NR>1 {{print $0, tool, bamtype, aligner, sample, coverage}}' "{output.smy}" > "{output.smy}.tmp" && \
            mv "{output.smy}.tmp" "{output.smy}"

        zcat {params.outdir}.roc.all.csv.gz | \
        awk -v tool="{wildcards.tool}" -v bamtype="{wildcards.bamtype}" -v aligner="{wildcards.aligner}" \
            -v sample="{wildcards.sample}" -v coverage="{wildcards.coverage}" \
        'BEGIN {{FS=OFS=","}}
            NR==1 {{print $0, "tool", "bamtype", "aligner", "sample", "coverage"}}
            NR>1 {{print $0, tool, bamtype, aligner, sample, coverage}}' | gzip > {output.roc}

        rm -f $TMP $TMP.tbi
        """


rule callable_bases_common:
    priority: 80
    input:
        "results/beds/common/{sample}_minimap2_{coverage}.bed"
    output:
        "results/beds/common/summary/{sample}_{coverage}.tsv"
    shell:
        r"""
        mkdir -p $(dirname {output})
        awk -v sample="{wildcards.sample}" -v mincov="{wildcards.coverage}" '
            BEGIN {{ OFS="\t"; callable=0 }}
            {{ callable += $3-$2 }}
            END {{ print "sample","status","min_coverage","metric","value";
                   print sample,"common",mincov,"callable_bases",callable }}' {input} > {output}
        """
