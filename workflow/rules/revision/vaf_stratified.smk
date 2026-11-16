# recall and precision stratified by pileup VAF, from hap.py output VCFs

VAF_TOOLS = ["Clair3-RNA", "longcallR-nn", "DeepVariant", "longcallR", "isoLASER", "GATK", "NanoTS"]


rule vaf_sites:
    input:
        vcfs=expand("results/happy/{tool}_origin_{{aligner}}_{{sample}}_{{coverage}}.vcf.gz", tool=VAF_TOOLS)
    output:
        "results/vaf/{sample}_{aligner}_{coverage}.sites.tsv"
    conda:
        "../../envs/bcftools.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        for v in {input.vcfs}; do
            bcftools query -f '%CHROM\\t%POS[\\t%BD\\t%BVT]\\n' $v \
                | awk -F'\\t' '(($3=="TP"||$3=="FN") && $4=="SNP") || ($5=="FP" && $6=="SNP") {{print $1"\\t"$2}}'
        done | sort -k1,1V -k2,2n -u > {output}
        """


rule vaf_pileup:
    input:
        bam="results/align/raw/origin/{aligner}_{sample}.aligned.bam",
        sites="results/vaf/{sample}_{aligner}_{coverage}.sites.tsv",
        ref=config["reference_genome"]
    output:
        "results/vaf/{sample}_{aligner}_{coverage}.vaf.tsv"
    conda:
        "../../envs/bcftools.yaml"
    params:
        threads=8,
        max_depth=2000,
        tmp="results/vaf/tmp_{sample}_{aligner}_{coverage}"
    resources:
        mem_mb=8000
    shell:
        """
        # VAF = non-reference reads / depth (MAPQ >= 20, no BAQ)
        mkdir -p {params.tmp}
        cut -f1 {input.sites} | sort -u \
            | xargs -P {params.threads} -I{{}} sh -c \
                'samtools mpileup -l {input.sites} -r {{}} -q 20 -Q 0 -B -d {params.max_depth} \
                    --no-output-ins --no-output-del --no-output-ends -f {input.ref} {input.bam} \
                    2> /dev/null > {params.tmp}/{{}}.pileup'
        cat {params.tmp}/*.pileup \
            | awk -F'\\t' 'BEGIN{{OFS="\\t"; print "chrom","pos","depth","nonref","vaf"}}
                {{s=toupper($5); gsub(/[$<>*]/,"",s); ref=gsub(/[.,]/,"",s); nr=length(s); d=ref+nr;
                  print $1,$2,d,nr,(d>0 ? nr/d : "NA")}}' > {output}
        rm -rf {params.tmp}
        """


rule vaf_labels:
    input:
        vcf="results/happy/{tool}_origin_{aligner}_{sample}_{coverage}.vcf.gz",
        vaf="results/vaf/{sample}_{aligner}_{coverage}.vaf.tsv"
    output:
        "results/vaf/labels/{tool}_{aligner}_{sample}_{coverage}.tsv"
    wildcard_constraints:
        tool="[^/_]+"
    conda:
        "../../envs/bcftools.yaml"
    shell:
        """
        mkdir -p $(dirname {output})
        bcftools query -f '%CHROM\\t%POS[\\t%BD\\t%BVT]\\n' {input.vcf} \
            | awk -F'\\t' -v OFS='\\t' -v tool="{wildcards.tool}" -v sample="{wildcards.sample}" \
                  -v aligner="{wildcards.aligner}" -v coverage="{wildcards.coverage}" '
                NR==FNR {{ if (FNR>1) {{ dp[$1":"$2]=$3; vf[$1":"$2]=$5 }}; next }}
                FNR==1 {{ print "tool","sample","aligner","coverage","chrom","pos","truth","query","depth","vaf" }}
                {{ t = (($3=="TP"||$3=="FN") && $4=="SNP") ? $3 : ".";
                   q = (($5=="TP"||$5=="FP") && $6=="SNP") ? $5 : ".";
                   if (t=="." && q==".") next;
                   k=$1":"$2;
                   print tool,sample,aligner,coverage,$1,$2,t,q,(k in dp ? dp[k] : "NA"),(k in vf ? vf[k] : "NA") }}' \
                {input.vaf} - > {output}
        """


ROBUST_TOOLS = ["Clair3-RNA", "DeepVariant", "isoLASER", "longcallR", "longcallR-nn", "NanoTS"]
ROBUST_SAMPLES = [f"{c}-{p}" for c in ["HG004", "HG005"]
                  for p in ["Baylor-IsoSeq", "MasSeq", "cDNAxR09", "cDNAxR10", "dRNA002", "dRNA004"]]

rule robustness_figure:
    input:
        script="code/revision/plt_robustness_GIAB.R",
        labels=expand("results/vaf/labels/{tool}_minimap2_{sample}_5.tsv", tool=ROBUST_TOOLS, sample=ROBUST_SAMPLES),
        roc=expand("results/happy/{tool}_origin_minimap2_{sample}_5.roc.csv.gz", tool=ROBUST_TOOLS, sample=ROBUST_SAMPLES)
    output:
        "plts/revision-robustness_GIAB.pdf"
    conda:
        "../../envs/r_env.yaml"
    log:
        "logs/plt_revision-robustness_GIAB.log"
    shell:
        "Rscript {input.script} {output} {input.labels} {input.roc} > {log} 2>&1"
