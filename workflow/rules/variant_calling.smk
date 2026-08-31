def get_clair3_rna_platform(wildcards):
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
    elif "r09" in sample and "cdna" in sample:
        return "ont_r9_guppy_cdna"
    else:
        raise ValueError(f"Cannot determine platform for sample '{wildcards.sample}' and aligner '{wildcards.aligner}'")

rule run_clair3_rna:
    priority: 95
    input:
        bam="results/align/raw/{bamtype}/{aligner}_{sample}.aligned.bam",
        ref=config["reference_genome"],
        img="/home/jiayiwang/tools/clair3-rna-latest.simg"
    output:
        dir=directory("/data/jiayiwang/variant-comparison/results/variant/Clair3-RNA/{bamtype}_{aligner}_{sample}"),
        file="results/variant/Clair3-RNA/{bamtype}_{aligner}_{sample}.vcf.gz"
    log: 
        "logs/run_clair3_rna_{bamtype},{aligner},{sample}.log"
    params:
        threads = 10,
        model = get_clair3_rna_platform
    resources:
        mem_mb=32000,
        tmpdir="/data/jiayiwang/variant-comparison/tmp",
    shell:
        """
        mkdir -p {resources.tmpdir}
        if [ ! -f {input[0]}.bai ]; then
            samtools index {input[0]};
        fi
        mkdir -p {output[0]}
        TMPDIR={resources.tmpdir} singularity exec \
            -B $(pwd)/{input[0]},$(pwd)/{input[1]} \
            -B /data/jiayiwang/variant-comparison:/data/jiayiwang/variant-comparison \
            {input[2]} \
            /bin/bash -c "source /opt/conda/bin/activate /opt/conda/envs/clair3_rna && \
            /opt/bin/run_clair3_rna \
            --bam_fn $(pwd)/{input[0]} \
            --ref_fn $(pwd)/{input[1]} \
            --threads {params.threads} \
            --platform {params.model} \
            --tag_variant_using_readiportal \
            --remove_intermediate_dir \
            --output_dir {output[0]} \
            --conda_prefix /opt/conda/envs/clair3_rna" \
            > {log} 2>&1

        if [ -z "$(zcat {output[0]}/output.vcf.gz | grep -v '^#' | head -1)" ]; then
            echo "ERROR: {output[0]}/output.vcf.gz has no variant records - run_clair3_rna likely failed silently, see {log}" >&2
            exit 1
        fi
        mv {output[0]}/output.vcf.gz {output[1]}
        """




rule run_deep_variant:
    priority: 95
    input: 
        bam=lambda wildcards: 
            f"results/align/preprocessed/{wildcards.bamtype}/{wildcards.aligner}_{wildcards.sample}.aligned.bam"
            if wildcards.sample in ISOSEQ_SAMPLE
            else f"results/align/raw/{wildcards.bamtype}/{wildcards.aligner}_{wildcards.sample}.aligned.bam",
        ref=config["reference_genome"],
        img="/home/jiayiwang/tools/deepvariant_1.9.0.sif"
    output:
        "results/variant/DeepVariant/{bamtype}_{aligner}_{sample}.vcf.gz",
    log:
        "logs/run_deep_variant_{bamtype},{aligner},{sample}.log"
    #threads: 10
    resources:
        mem_mb=32000,
        tmpdir="/data/jiayiwang/variant-comparison/tmp",
    params:
        threads = 10,
        model = lambda wildcards: "MASSEQ" if "IsoSeq" in wildcards.sample or "MasSeq" in wildcards.sample else "ONT_R104",
        extra_args = lambda wildcards: "" if "IsoSeq" in wildcards.sample or "MasSeq" in wildcards.sample else '--make_examples_extra_args="split_skip_reads=true"'
    shell:
        """
        mkdir -p {resources.tmpdir}
        if [ ! -f {input.bam}.bai ]; then
        samtools index {input.bam};
        fi
        export OMP_NUM_THREADS={threads}
        export TF_NUM_INTEROP_THREADS=1
        export TF_NUM_INTRAOP_THREADS={threads}
        export TMPDIR={resources.tmpdir}
        singularity exec --bind /usr/lib/locale/ --bind $(pwd -P):$(pwd -P) --pwd $(pwd -P) {input.img} /opt/deepvariant/bin/run_deepvariant\
            --model_type {params.model} \
            --ref {input.ref} \
            --reads {input.bam} \
            --output_vcf {output} \
            --num_shards {params.threads} \
            --intermediate_results_dir /home/jiayiwang/tmp/{wildcards.sample} \
        > {log} 2>&1
        """

def choose_pcr_model(wildcards):
    s = wildcards.sample
    if "dRNA" in s:
        return "NONE"
    elif "cDNA" in s or "IsoSeq" in s or "MasSeq" in s:
        return "CONSERVATIVE"
    else:
        return "NONE"

rule run_GATK:
    priority: 95
    input:
        bam = "results/align/raw/{bamtype}/{aligner}_{sample}.aligned.bam",
        img = "/home/jiayiwang/tools/gatk-latest.simg"
    output:
        vcf = "results/variant/GATK/{bamtype}_{aligner}_{sample}.vcf.gz"
    conda:
        "../envs/bcftools.yaml"
    params:
        ref = config["reference_genome"],
        threads = 10,
        pcr_model = choose_pcr_model
    resources:
        mem_mb=10000,
        tmpdir="/data/jiayiwang/variant-comparison/tmp",
    shell:
        """
        mkdir -p {resources.tmpdir}
        tmpbam=$(mktemp --tmpdir={resources.tmpdir} --suffix=.bam)
        tmpsrt=$(mktemp --tmpdir={resources.tmpdir} --suffix=.bam)
        samtools addreplacerg \
            -r '@RG\\tID:1\\tSM:{wildcards.sample}\\tPL:ONT' \
            -o "$tmpbam" {input.bam}
        
        samtools sort -@ {params.threads} -o "$tmpsrt" "$tmpbam"
        rm "$tmpbam"
        samtools index "$tmpsrt"

        singularity exec --bind $(pwd -P):$(pwd -P) --pwd $(pwd -P) --bind /home/jiayiwang/miniconda3/envs/snakemake/bin/python:/usr/bin/python \
        {input.img} /gatk/gatk --java-options "-Xmx4G" HaplotypeCaller \
            -R {params.ref} \
            -I "$tmpsrt" \
            -O {output.vcf} \
            --dont-use-soft-clipped-bases true \
            --pcr-indel-model {params.pcr_model} \
            --native-pair-hmm-threads {params.threads}

        rm "$tmpsrt" "$tmpsrt.bai"
        if zgrep -q 'UnnamedSample' "{output.vcf}"; then
            echo "Fixing malformed VCF header..."
            bcftools view -s $(bcftools query -l "{output.vcf}" | head -n1) -Oz -o "{output.vcf}.tmp" "{output.vcf}" && \
            mv "{output.vcf}.tmp" "{output.vcf}" && \
            bcftools index -t -f "{output.vcf}"
        fi
        """

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

CHR = [f"chr{i}" for i in range(1, 23)]

rule run_longcallR_chr:
    priority: 95
    input:
        "results/align/raw/{bamtype}/{aligner}_{sample}.aligned.bam"
    output:
        vcf="results/variant/longcallR/{bamtype}_{aligner}_{sample}/vcf/{chr}.vcf"
    log:
        "logs/run_longcallR_{bamtype}_{aligner}_{sample}_{chr}.log"
    params:
        threads = 10,
        platform = get_longcallR_platform,
        ref = config["reference_genome"]
    resources:
        mem_mb=16000,
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        if [ ! -f {input}.bai ] || [ {input} -nt {input}.bai ]; then
            samtools index {input};
        fi

        mkdir -p $(dirname {output.vcf})

        ~/tools/longcallR/longcallR \
            --bam-path {input} \
            --ref-path {params.ref} \
            --output $(dirname {output.vcf})/{wildcards.chr} \
            --preset {params.platform} \
            --contigs {wildcards.chr} \
            -t {params.threads} \
            --no-bam-output \
            --min-depth 0 \
            --max-depth 50000 \
            > {log} 2>&1
        """

rule run_longcallR:
    priority: 95
    input:
        vcfs=expand("results/variant/longcallR/{{bamtype}}_{{aligner}}_{{sample}}/vcf/{chr}.vcf", chr=CHR)
    output:
        file="results/variant/longcallR/{bamtype}_{aligner}_{sample}.vcf.gz"
    log:
        "logs/run_longcallR_merge_{bamtype}_{aligner}_{sample}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        (
          bcftools concat {input.vcfs} 2>> {log} \
          | bcftools sort -Oz -o {output.file} 2>> {log}
          tabix -p vcf {output.file} 2>> {log}
        ) >> {log} 2>&1
        """


def pick_longcallr_nn_model(sample, model_dir="/home/jiayiwang/tools/longcallR-nn/models"):
    s = str(sample).lower()
    if "cdna" in s:
        return (f"{model_dir}/ont_cdna_config.yaml", f"{model_dir}/ont_cdna_model.chkpt")
    if "drna" in s:
        return (f"{model_dir}/ont_drna_config.yaml", f"{model_dir}/ont_drna_model.chkpt")
    if "isoseq" in s:
        return (f"{model_dir}/pb_isoseq_config.yaml", f"{model_dir}/pb_isoseq_model.chkpt")
    if "masseq" in s:
        return (f"{model_dir}/pb_masseq_config.yaml", f"{model_dir}/pb_masseq_model.chkpt")
    raise ValueError(f"sample '{sample}' does not contain cDNA/dRNA/IsoSeq/MasSeq")

rule longcallR_dp_predict:
    priority: 95
    input:
        bam="results/align/raw/{bamtype}/{aligner}_{sample}.aligned.bam"
    output:
        dir=directory("results/variant/longcallR-nn/{bamtype}_{aligner}_{sample}/data/{chr}")
    log:
        "logs/longcallR_dp_predict_{bamtype}_{aligner}_{sample}_{chr}.log"
    params:
        ref=config["reference_genome"],
        thread=10,
        min_baseq=lambda wc: "" if ("MasSeq" in wc.sample or "IsoSeq" in wc.sample) else "--min-baseq 10"
    resources:
        mem_mb=12000,
    shell:
        """
        if [ ! -f {input.bam}.bai ] || [ {input.bam} -nt {input.bam}.bai ]; then
            samtools index {input.bam};
        fi
        mkdir -p {output.dir}
        /home/jiayiwang/tools/longcallR-nn/longcallR_dp/target/release/longcallR_dp --mode predict \
          --bam-path {input.bam} \
          --ref-path {params.ref} \
          --threads {params.thread} \
          --contigs {wildcards.chr} \
          --output {output.dir} \
          {params.min_baseq} \
          >> {log} 2>&1
        """


rule longcallR_nn_call:
    priority: 95
    input:
        data="results/variant/longcallR-nn/{bamtype}_{aligner}_{sample}/data/{chr}",
        bam="results/align/raw/{bamtype}/{aligner}_{sample}.aligned.bam"
    output:
        vcf="results/variant/longcallR-nn/{bamtype}_{aligner}_{sample}/vcf/{chr}.vcf"
    log:
        "logs/longcallR_nn_call_{bamtype}_{aligner}_{sample}_{chr}.log"
    params:
        ref=config["reference_genome"],
        no_cuda=True,
        cfg=lambda wc: pick_longcallr_nn_model(wc.sample)[0],
        ckpt=lambda wc: pick_longcallr_nn_model(wc.sample)[1]
    resources:
        mem_mb=12000,
    shell:
        """
        NO_CUDA=""
        if [ "{params.no_cuda}" = "True" ]; then NO_CUDA="--no_cuda"; fi

        mkdir -p $(dirname {output.vcf})

        /home/jiayiwang/miniconda3/envs/longcallRenv/bin/longcallR_nn call \
          -config "{params.cfg}" \
          -model "{params.ckpt}" \
          -data "{input.data}" \
          -ref {params.ref} \
          -output {output.vcf} \
          $NO_CUDA \
          -max_depth 200 \
          -batch_size 256 \
          >> {log} 2>&1
        """

rule longcallR_nn_merge:
    priority: 95
    input:
        vcfs=expand("results/variant/longcallR-nn/{{bamtype}}_{{aligner}}_{{sample}}/vcf/{chr}.vcf", chr=CHR)
    output:
        vcfgz="results/variant/longcallR-nn/{bamtype}_{aligner}_{sample}.vcf.gz",
        tbi="results/variant/longcallR-nn/{bamtype}_{aligner}_{sample}.vcf.gz.tbi"
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/longcallR_nn_merge_{bamtype}_{aligner}_{sample}.log"
    shell:
        """
        (
          bcftools concat {input.vcfs} 2>> {log} \
          | bcftools sort -Oz -o {output.vcfgz} 2>> {log}
          tabix -p vcf {output.vcfgz} 2>> {log}
        ) >> {log} 2>&1
        """

rule extract_transcriptome:
    priority: 95
    input:
        gtf=config["gtfgz"],
        ref=config["reference_genome"]
    output:
        fa="data/isolaser/transcriptome.fa",
        db=directory("data/isolaser/transcriptome.db"),
        gtf="data/isolaser/filtered.gtf.gz",
        tbi="data/isolaser/filtered.gtf.gz.tbi"
    conda:
        "../envs/bcftools.yaml"
    resources:
        mem_mb=6000,
    log:
        "logs/extract_transcriptome.log"
    shell:
        r"""
        zcat {input.gtf} | awk '$1 !~ /\./' | bgzip > {output.gtf}
        tabix -p gff {output.gtf}
        isolaser_convert_gtf_to_fasta -g {output.gtf} -f {input.ref} -o {output.fa} 2>> {log}
        mkdir -p {output.db}
        isolaser_extract_exon_parts -g {output.gtf} -o {output.db} 2>> {log}
        """


rule isolaser_realign_annotate:
    priority: 95
    input:
        bam=lambda wildcards: 
            f"results/align/preprocessed/{wildcards.bamtype}/{wildcards.aligner}_{wildcards.sample}.aligned.bam"
            if wildcards.sample in ISOSEQ_SAMPLE
            else f"results/align/raw/{wildcards.bamtype}/{wildcards.aligner}_{wildcards.sample}.aligned.bam",
        tx="data/isolaser/transcriptome.fa",
        gtf="data/isolaser/filtered.gtf.gz",
    output:
        bam="results/variant/isoLASER/realign/{bamtype}/{aligner}_{sample}.anno.bam"
    log:
        "logs/isolaser_annotate_{bamtype},{aligner},{sample}.log"
    params:
        sam="tmp/{bamtype}/{aligner}_{sample}.aligned.sam",
        unsorted="results/variant/isoLASER/realign/{bamtype}/{aligner}_{sample}.anno.unsorted.bam",
        threads=25,
    resources:
        mem_mb=40000,
        tmpdir="/data/jiayiwang/variant-comparison/tmp",
    conda:
        "../envs/isolaser.yaml"
    shell:
        """
        mkdir -p $(dirname {params.sam})
        mkdir -p $(dirname {output.bam})

        samtools bam2fq -@ {params.threads} {input.bam} \
            | minimap2 -t {params.threads} -ax splice:hq -uf --MD {input.tx} - \
            > {params.sam} 2>> {log}

        isolaser_annotate \
            -b {input.bam} \
            -t {params.sam} \
            -g {input.gtf} \
            -o {params.unsorted} 2>> {log}

        samtools sort -@ {params.threads} -m 4G \
            -T {resources.tmpdir}/{wildcards.sample}_sort \
            -o {output.bam} {params.unsorted} 2>> {log}

        samtools index {output.bam} 2>> {log}

        rm -f {params.sam} {params.unsorted}
        """

def get_isoLASER_platform(wildcards):
    sample = wildcards.sample.lower()
    for key, platform in {
        "masseq": "PacBio",
        "isoseq": "PacBio",
        "cdna": "Nanopore",
        "drna": "Nanopore"
    }.items():
        if key in sample:
            return platform
    raise ValueError(f"Unrecognized platform for sample '{wildcards.sample}'")

rule run_isoLASER:
    priority: 95
    input:
        bam=rules.isolaser_realign_annotate.output.bam,
        gtf="data/isolaser/filtered.gtf.gz",
        ref=config["reference_genome"],
        db="data/isolaser/transcriptome.db"
    output:
        vcfgz="results/variant/isoLASER/{bamtype}_{aligner}_{sample}.vcf.gz",
        tbi="results/variant/isoLASER/{bamtype}_{aligner}_{sample}.vcf.gz.tbi"
    log:
        "logs/run_isoLASER_{bamtype},{aligner},{sample}.log"
    conda:
        "../envs/isolaser.yaml"
    resources:
        mem_mb=40000
    params:
        threads=25,
        platform = get_isoLASER_platform,
        prefix="results/variant/isoLASER/{bamtype}_{aligner}_{sample}/calls"
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        ref="{input.ref}"
        clean_ref="${{ref%.fa}}.clean.fa"
        if [ ! -f "$clean_ref" ]; then
            awk '/^>/{{print; next}} {{gsub(/[^ACGTNacgtn]/, "N"); print}}' {input.ref} > "$clean_ref"
        fi
        isolaser \
            -b {input.bam} \
            -o {params.prefix} \
            -t {input.db} \
            -f "$clean_ref" \
            --DP=0 \
            -n {params.threads} \
            --platform={params.platform} \
            --include-continuous \
            > {log} 2>&1
        bcftools view -v snps,indels {params.prefix}.gvcf \
            | bcftools norm -m -any \
            | bcftools view -e 'GT="0/0"' -Oz \
            -o {output.vcfgz}
        tabix -p vcf {output.vcfgz}
        """


rule isolaser_filter:
    input:
        vcf = "results/variant/isoLASER/{bamtype}_{aligner}_{sample}.vcf.gz"
    output:
        vcf = "results/variant/isoLASER_filtered/{bamtype}_{aligner}_{sample}.vcf.gz"
    log:
        "logs/isolaser_filter_{bamtype}_{aligner}_{sample}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        mkdir -p $(dirname {output.vcf})

        # Author-provided evaluation criteria: allelic ratio (AB) > 25%, GQ > 25,
        # applied to all variant types. For multi-ALT sites, keep only the top-2 alleles by
        # quality after splitting to biallelic records.
        bcftools view -e 'ALT="<NON_REF>"' -f PASS {input.vcf} \
            | bcftools norm -m -any \
            | bcftools view -e 'ALT="<NON_REF>"' \
            | bcftools view -i 'FMT/GQ>25 && FMT/AB>0.25' \
            | awk -F'\t' 'BEGIN{{OFS="\t"; n=0}}
                /^#/ {{print; next}}
                {{
                    key=$1"\t"$2
                    if (key!=prevkey && prevkey!="") flush()
                    buf[n]=$0; qual[n]=$6+0; n++
                    prevkey=key
                }}
                END{{ if(n>0) flush() }}
                function flush(   i,j,tmp){{
                    for(i=0;i<n;i++) for(j=i+1;j<n;j++) if(qual[j]>qual[i]){{tmp=buf[i];buf[i]=buf[j];buf[j]=tmp; tmp=qual[i];qual[i]=qual[j];qual[j]=tmp}}
                    for(i=0;i<n && i<2;i++) print buf[i]
                    n=0; delete buf; delete qual
                }}' \
            | bcftools norm --rm-dup both -Oz -o {output.vcf} \
            2> {log}

        bcftools index -t -f {output.vcf} 2>> {log}
        """

