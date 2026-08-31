rule generate_bed:
    priority: 50
    input: 
        bam="results/align/{status}/origin/{aligner}_{sample}.aligned.bam",
        gtf=config["gtf"]
    output:
        "results/beds/{status}/{sample}_{aligner}_{coverage}.bed"
    conda:
        "../envs/bcftools.yaml"
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
        singularity exec \
            --bind $(pwd -P):$(pwd -P) \
            --pwd $(pwd -P) \
            {params.img} mosdepth \
            --threads {params.threads} \
            results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage} \
            {input.bam}

        gzip -fdc results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.per-base.bed.gz \
            | awk -v min_cov={params.min_cov} '$4 >= min_cov' \
            | bedtools merge -d 1 -c 4 -o mean -i - \
            > results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.coverage.bed

        bedtools intersect \
            -a results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.coverage.bed \
            -b <(grep -v "^#" {input.gtf} | awk '$3 == "exon"') \
            > results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.tmp.bed

        if [ -f {params.confidence} ]; then
            bedtools intersect \
                -a results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.tmp.bed \
                -b {params.confidence} > {output}
            rm results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.tmp.bed
        else
            mv results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}.tmp.bed {output}
        fi
        rm -rf results/beds/{wildcards.status}/{wildcards.sample}/{wildcards.aligner}_{wildcards.coverage}
        """

rule callable_bases:
    priority: 80
    input:
        "results/beds/raw/{sample}_minimap2_{coverage}.bed"
    output:
        "results/beds/raw/summary/{sample}_{coverage}.tsv"
    shell:
        r"""
        mkdir -p $(dirname {output})

        awk -v sample="{wildcards.sample}" \
            -v status="raw" \
            -v mincov="{wildcards.coverage}" '
        BEGIN{{
            OFS="\t";
            callable=0;
            weighted_cov=0;
        }}
        {{
            len=$3-$2;
            callable += len;
            weighted_cov += len*$4;
        }}
        END{{
            mean_cov = (callable>0 ? weighted_cov/callable : "NA");

            print "sample","status","min_coverage","metric","value";
            print sample,status,mincov,"callable_bases",callable;
            print sample,status,mincov,"mean_callable_coverage",mean_cov;
            print sample,status,mincov,"covered_bases",weighted_cov;
        }}' {input} > {output}
        """

# rule callable_bases_context:
#     priority: 80
#     input:
#         bed="results/beds/raw/{sample}_minimap2_{coverage}.bed",
#         region=lambda wildcards: STRAT_BEDS[wildcards.subset]
#     output:
#         "results/beds/context/{sample}_{coverage}_{subset}.tsv"
#     wildcard_constraints:
#         coverage=r"\d+",
#         subset="|".join(CONTEXT_SUBSETS)
#     conda:
#         "../envs/bcftools.yaml"
#     shell:
#         r"""
#         mkdir -p $(dirname {output})

#         bedtools intersect -a {input.bed} -b {input.region} \
#             | awk -v sample="{wildcards.sample}" \
#                   -v mincov="{wildcards.coverage}" \
#                   -v subset="{wildcards.subset}" '
#         BEGIN{{OFS="\t"; callable=0}}
#         {{callable += $3-$2}}
#         END{{
#             print "sample","coverage","subset","metric","value";
#             print sample,mincov,subset,"callable_bases",callable;
#         }}' > {output}
#         """

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


def get_variant_vcf(wildcards):
    # isoLASER's GQ/AB/multi-ALT filtering lives in its own rule (isolaser_filter,
    # variant_calling.smk) rather than inline here; every other tool is used as-is.
    if wildcards.tool == "isoLASER":
        return f"results/variant/isoLASER_filtered/{wildcards.bamtype}_{wildcards.aligner}_{wildcards.sample}.vcf.gz"
    return f"results/variant/{wildcards.tool}/{wildcards.bamtype}_{wildcards.aligner}_{wildcards.sample}.vcf.gz"


rule happy_benchmark:
    priority: 0
    retries: 2
    input:
        vcf=get_variant_vcf,
        bed=get_bed,
        confidence = lambda wildcards: f"data/truth/beds/{wildcards.sample.split('-')[0].rstrip('abc')}.bed"
    output:
        smy="results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.summary.csv",
        roc="results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.roc.csv.gz",
        vcf="results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    params:
        truth = lambda wildcards: f"data/truth/{wildcards.sample.split('-')[0].rstrip('abc')}/truth.vcf.gz",
        ref = config["reference_genome"],
        outdir = "results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}",
        threads = 10 
    log:
        "logs/happy_{tool}_{bamtype}_{aligner}_{sample}_{coverage}.log"
    shell:
        """
        export HGREF={params.ref}

        JOBTAG={wildcards.tool}_{wildcards.bamtype}_{wildcards.aligner}_{wildcards.sample}_{wildcards.coverage}
        SRC=tmp/${{JOBTAG}}.src.vcf.gz
        cp "{input.vcf}" "$SRC"

        if ! bcftools index -t -f "$SRC" >/dev/null 2>&1; then
            echo "VCF is not sorted or not indexed properly. Sorting and indexing..."
            mv "$SRC" "$SRC.unsorted"
            bcftools sort -Oz -o "$SRC" "$SRC.unsorted"
            bcftools index -t -f "$SRC"
            rm "$SRC.unsorted"
        fi

        TMP=tmp/${{JOBTAG}}.vcf.gz

        if [ "{wildcards.tool}" = "GATK" ]; then
            echo "GATK detected → using QUAL>=30" >&2
            bcftools view -i 'QUAL>=30' -Oz -o $TMP "$SRC"

        else
            echo "Non-GATK → using FILTER=PASS" >&2
            bcftools view -f PASS -Oz -o $TMP "$SRC"
        fi

        bcftools index -t -f $TMP
        rm -f "$SRC" "$SRC.tbi"

        singularity exec \
            --bind $(pwd -P):$(pwd -P) \
            --pwd $(pwd -P) \
            /home/jiayiwang/tools/hap.py.simg /opt/hap.py/bin/hap.py \
            {params.truth} \
            $TMP \
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

        rm $TMP
        rm -f $TMP.tbi
        """




rule stratified_benchmark:
    input:
        vcf=rules.happy_benchmark.output.vcf,
        stratification="data/regions/GRCh38-all-stratifications.tsv"
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
            --bind $(pwd):$(pwd) \
            --pwd $(pwd) \
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
    conda:
        "../envs/snpeff.yaml"
    shell:
        """
        snpEff -Xmx8g GRCh38.86 {input.vcf} | bgzip -c > {output}
        """

rule redi_bed:
    input:
        "data/rna_editing/TABLE1_hg38_v3.txt.gz"
    output:
        "data/rna_editing/REDIportal.bed.gz"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        # id column encodes the site's own Ref>Ed alleles (genome-forward-strand,
        # per REDIportal's Ref/Ed columns) in the same CHROM_POS_REF_ALT format used
        # for VCF calls elsewhere, so downstream rules can require an exact
        # substitution match rather than just positional overlap (REDIportal is
        # exclusively A>G / T>C, but a caller's FP at that position may not be).
        zcat {input} | tail -n +2 | awk -F'\t' '{{print $2"\t"($3-1)"\t"$3"\t"$2"_"$3"_"$4"_"$5}}' \
            | sort -k1,1 -k2,2n \
            | uniq \
            | bgzip > {output}
        """

rule callable_redi_sites:
    input:
        bed="results/beds/raw/{sample}_minimap2_{coverage}.bed",
        redi=rules.redi_bed.output
    output:
        "results/rna_editing/callable/{sample}_{coverage}.tsv"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        mkdir -p $(dirname {output})

        sort -k1,1 -k2,2n {input.bed} > {output}.sorted.bed

        n=$(bedtools intersect -u -a {input.redi} -b {output}.sorted.bed | wc -l)

        printf 'sample\tcoverage\tcallable_redi_sites\n%s\t%s\t%s\n' \
            "{wildcards.sample}" "{wildcards.coverage}" "$n" > {output}

        rm -f {output}.sorted.bed
        """

rule rna_editing_overlap:
    input:
        vcf="results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.vcf.gz",
        orig_vcf=get_variant_vcf,
        redi=rules.redi_bed.output
    output:
        "results/rna_editing/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.tsv"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        mkdir -p $(dirname {output})

        bcftools view -f PASS {input.orig_vcf} \
            | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
            | awk -F'\t' 'BEGIN{{OFS="\t"}} {{print $1,$2-1,$2-1+length($3),$1"_"$2"_"$3"_"$4}}' \
            | sort -k1,1 -k2,2n > {output}.pass.bed

        bcftools query -s QUERY -f '%CHROM\t%POS\t%REF\t%ALT\t[%BD]\t[%BVT]\n' {input.vcf} \
            | awk -F'\t' '$5=="FP"{{print $1"\t"($2-1)"\t"($2-1+length($3))"\t"$1"_"$2"_"$3"_"$4"\t"$6}}' \
            | sort -k1,1 -k2,2n > {output}.fp.bed

        awk -F'\t' 'NR==FNR {{pass[$4]=1; next}} ($4 in pass)' \
            {output}.pass.bed {output}.fp.bed > {output}.fp_pass.bed

        # Loose: genomic position overlap only -- any FP landing on a cataloged
        # editing coordinate counts, regardless of its own substitution.
        bedtools intersect -u -a {output}.fp_pass.bed -b {input.redi} > {output}.fp_redi_loose.bed

        # Strict: exact CHROM_POS_REF_ALT match against REDIportal's Ref>Ed alleles.
        # REDIportal's Ref/Ed are already given in genome-forward-strand coordinates
        # (strand-corrected per site: e.g. a minus-strand gene's A-to-I edit is
        # recorded as T>C, not A>G), so matching a FP's own REF>ALT against the
        # specific entry at that position enforces both the A-to-I edit signature
        # and strand-consistency simultaneously.
        zcat {input.redi} | cut -f4 | sort -u > {output}.redi_ids.txt
        awk -F'\t' 'NR==FNR {{redi[$1]=1; next}} ($4 in redi)' \
            {output}.redi_ids.txt {output}.fp_pass.bed > {output}.fp_redi_strict.bed

        # Per-file type (SNP/INDEL) counts, always emitting both keys even if zero,
        # so joining below is safe regardless of which files are empty.
        for f in {output}.fp_pass.bed {output}.fp_redi_loose.bed {output}.fp_redi_strict.bed; do
            awk -F'\t' 'BEGIN{{c["SNP"]=0;c["INDEL"]=0}} {{c[$5]++}} END{{for (t in c) print t"\t"c[t]}}' "$f" > "$f.counts"
        done

        awk -v tool="{wildcards.tool}" -v bamtype="{wildcards.bamtype}" \
            -v aligner="{wildcards.aligner}" -v sample="{wildcards.sample}" \
            -v coverage="{wildcards.coverage}" -F'\t' '
        BEGIN {{OFS="\t"; print "tool","bamtype","aligner","sample","coverage","Type","FP","FP_in_REDIportal_loose","FP_in_REDIportal_strict"}}
        FNR==1 {{ file++ }}
        file==1 {{ fp[$1]=$2; next }}
        file==2 {{ loose[$1]=$2; next }}
        file==3 {{ strict[$1]=$2; next }}
        END {{
            for (t in fp) print tool,bamtype,aligner,sample,coverage,t,fp[t],loose[t],strict[t]
        }}' {output}.fp_pass.bed.counts {output}.fp_redi_loose.bed.counts {output}.fp_redi_strict.bed.counts > {output}

        rm -f {output}.pass.bed {output}.fp.bed {output}.fp_pass.bed \
            {output}.fp_redi_loose.bed {output}.fp_redi_strict.bed {output}.redi_ids.txt \
            {output}.fp_pass.bed.counts {output}.fp_redi_loose.bed.counts {output}.fp_redi_strict.bed.counts
        """


rule happy_chr20:
    input:
        vcf="results/happy/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.vcf.gz"
    output:
        "results/happy_chr20/{tool}_{bamtype}_{aligner}_{sample}_{coverage}.tsv"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        mkdir -p $(dirname {output})

        bcftools query -r chr20 -s TRUTH,QUERY -f '%CHROM\t%POS\t[%BD,]\t[%BVT,]\n' {input.vcf} \
            | awk -v tool="{wildcards.tool}" -v bamtype="{wildcards.bamtype}" \
                  -v aligner="{wildcards.aligner}" -v sample="{wildcards.sample}" \
                  -v coverage="{wildcards.coverage}" -F'\t' '
        BEGIN {{OFS="\t"; print "tool","bamtype","aligner","sample","coverage","Type","TP","FP","FN"}}
        {{
            split($3, bd, ",");
            split($4, bvt, ",");
            truth_bd = bd[1]; query_bd = bd[2];
            truth_bvt = bvt[1]; query_bvt = bvt[2];
            if (truth_bd == "TP") {{ tp[truth_bvt]++; types[truth_bvt] = 1 }}
            if (truth_bd == "FN") {{ fn[truth_bvt]++; types[truth_bvt] = 1 }}
            if (query_bd == "FP") {{ fp[query_bvt]++; types[query_bvt] = 1 }}
        }}
        END {{
            for (t in types) print tool,bamtype,aligner,sample,coverage,t,tp[t]+0,fp[t]+0,fn[t]+0
        }}' > {output}
        """