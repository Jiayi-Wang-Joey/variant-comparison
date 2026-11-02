# hap.py after removing calls matching REDIportal sites


rule redi_vcf:
    input:
        table="data/rna_editing/TABLE1_hg38_v3.txt.gz",
        fai=config["reference_genome"] + ".fai"
    output:
        vcf="data/rna_editing/REDIportal.sites.vcf.gz",
        tbi="data/rna_editing/REDIportal.sites.vcf.gz.tbi"
    conda:
        "../../envs/bcftools.yaml"
    shell:
        """
        mkdir -p /data/jiayiwang/tmp
        (
            echo '##fileformat=VCFv4.2'
            awk '{{print "##contig=<ID="$1",length="$2">"}}' {input.fai}
            printf '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n'
            zcat {input.table} | tail -n +2 \
                | awk -F'\\t' -v OFS='\\t' 'NR==FNR {{ok[$1]=1; next}} ($2 in ok) {{print $2,$3,$1,$4,$5,".",".","."}}' \
                    {input.fai} - \
                | sort -k1,1 -k2,2n -u
        ) | bcftools sort -T /data/jiayiwang/tmp/redi_sort.XXXXXX -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule happy_redimasked:
    priority: 0
    retries: 2
    input:
        vcf=get_variant_vcf,
        bed=get_bed,
        redi="data/rna_editing/REDIportal.sites.vcf.gz",
        confidence=lambda wc: f"data/truth/beds/{wc.sample.split('-')[0].rstrip('abc')}.bed"
    output:
        smy="results/happy_redimasked/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.summary.csv",
        roc="results/happy_redimasked/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.roc.csv.gz",
        vcf="results/happy_redimasked/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.vcf.gz",
        removed="results/happy_redimasked/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.n_removed.tsv"
    wildcard_constraints:
        tool="[^/_]+"
    conda:
        "../../envs/bcftools.yaml"
    params:
        truth=lambda wc: f"data/truth/{wc.sample.split('-')[0].rstrip('abc')}/truth.vcf.gz",
        ref=config["reference_genome"],
        outdir="results/happy_redimasked/{tool}_{bamtype}_{aligner}_{sample}_{coverage}",
        threads=10
    log:
        "logs/happy_redimasked_{tool}_{bamtype}_{aligner}_{sample}_{coverage}.log"
    shell:
        """
        export HGREF={params.ref}
        mkdir -p results/happy_redimasked tmp /data/jiayiwang/tmp

        JOBTAG=redimasked_{wildcards.tool}_{wildcards.bamtype}_{wildcards.aligner}_{wildcards.sample}_{wildcards.coverage}
        CALLS=tmp/${{JOBTAG}}.calls.vcf.gz
        TMP=tmp/${{JOBTAG}}.vcf.gz

        # GATK QUAL>=30, others PASS; split multi-allelic
        if [ "{wildcards.tool}" = "GATK" ]; then
            bcftools view -i 'QUAL>=30' {input.vcf} | bcftools sort -T /data/jiayiwang/tmp/${{JOBTAG}}.XXXXXX | bcftools norm -m -any -Oz -o $CALLS
        else
            bcftools view -f PASS {input.vcf} | bcftools sort -T /data/jiayiwang/tmp/${{JOBTAG}}.XXXXXX | bcftools norm -m -any -Oz -o $CALLS
        fi
        bcftools index -t -f $CALLS

        bcftools isec -C -c none -w1 $CALLS {input.redi} -Oz -o $TMP
        bcftools index -t -f $TMP
        printf 'tool\\tsample\\tcoverage\\tcalls\\tremoved_redi\\n%s\\t%s\\t%s\\t%s\\t%s\\n' \
            "{wildcards.tool}" "{wildcards.sample}" "{wildcards.coverage}" \
            "$(bcftools view -H $CALLS | wc -l)" \
            "$(( $(bcftools view -H $CALLS | wc -l) - $(bcftools view -H $TMP | wc -l) ))" > {output.removed}

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

        rm -f $CALLS $CALLS.tbi $TMP $TMP.tbi
        """


rule redimasked_combined:
    # single table: too many summaries to pass as arguments
    input:
        masked=lambda wc: happy_redimasked,
        orig=lambda wc: [f.replace("happy_redimasked", "happy") for f in happy_redimasked]
    output:
        "results/happy_redimasked/combined_SNP.csv"
    run:
        import csv
        keep = ["tool", "sample", "coverage", "TRUTH.TP", "TRUTH.FN", "QUERY.FP"]
        with open(output[0], "w", newline="") as out:
            w = csv.writer(out)
            w.writerow(["set"] + keep)
            for label, files in [("orig", input.orig), ("masked", input.masked)]:
                for f in files:
                    with open(f) as fh:
                        for r in csv.DictReader(fh):
                            if r["Type"] == "SNP" and r["Filter"] == "PASS":
                                w.writerow([label] + [r[k] for k in keep])
