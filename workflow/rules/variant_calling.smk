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
        dir=directory("/home/jiayiwang/variant-comparison/results/variant/Clair3-RNA/{bamtype}_{aligner}_{sample}"),
        file="results/variant/Clair3-RNA/{bamtype}_{aligner}_{sample}.vcf.gz"
    log: 
        "logs/run_clair3_rna_{bamtype},{aligner},{sample}.log"
    params:
        threads = 10,
        model = get_clair3_rna_platform
    shell:
        """
        if [ ! -f {input[0]}.bai ]; then
            samtools index {input[0]};
        fi
        singularity exec -B /home/jiayiwang/variant-comparison/{input[0]},{input[1]} {input[2]} \
            /bin/bash -c "source /opt/conda/bin/activate /opt/conda/envs/clair3_rna && \
            /opt/bin/run_clair3_rna \
            --bam_fn /home/jiayiwang/variant-comparison/{input[0]} \
            --ref_fn {input[1]} \
            --threads {params.threads} \
            --platform {params.model} \
            --tag_variant_using_readiportal \
            --remove_intermediate_dir \
            --output_dir {output[0]} \
            --conda_prefix /opt/conda/envs/clair3_rna" \
            > {log} 2>&1
        mv {output[0]}/output.vcf.gz {output[1]}
        """

# rule run_clair3_rna_phasing:
#     priority: 95
#     input:
#         bam="results/align/raw/{bamtype}/{aligner}_{sample}.aligned.bam",
#         ref=config["reference_genome"],
#         img="/home/jiayiwang/tools/clair3-rna-latest.simg"
#     output:
#         dir=directory("/home/jiayiwang/variant-comparison/results/variant/Clair3-RNA/{bamtype}_{aligner}_{sample}"),
#         file="results/variant/Clair3-RNA/{bamtype}_{aligner}_{sample}.vcf.gz"
#     log: 
#         "logs/run_clair3_rna_{bamtype},{aligner},{sample}.log"
#     params:
#         threads = 10,
#         model = get_clair3_rna_platform
#     shell:
#         """
#         if [ ! -f {input[0]}.bai ]; then
#             samtools index {input[0]};
#         fi
#         singularity exec -B /home/jiayiwang/variant-comparison/{input[0]},{input[1]} {input[2]} \
#             /bin/bash -c "source /opt/conda/bin/activate /opt/conda/envs/clair3_rna && \
#             /opt/bin/run_clair3_rna \
#             --bam_fn /home/jiayiwang/variant-comparison/{input[0]} \
#             --ref_fn {input[1]} \
#             --threads {params.threads} \
#             --platform {params.model} \
#             --tag_variant_using_readiportal \
#             --remove_intermediate_dir \
#             --output_dir {output[0]} \
#             --conda_prefix /opt/conda/envs/clair3_rna" \
#             > {log} 2>&1
#         mv {output[0]}/output.vcf.gz {output[1]}
#         """



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
    threads: 10
    params: 
        threads = 10,
        model = lambda wildcards: "MASSEQ" if "IsoSeq" in wildcards.sample or "MasSeq" in wildcards.sample else "ONT_R104",
        extra_args = lambda wildcards: "" if "IsoSeq" in wildcards.sample or "MasSeq" in wildcards.sample else '--make_examples_extra_args="split_skip_reads=true"'
    shell:
        """
        if [ ! -f {input.bam}.bai ]; then
        samtools index {input.bam};
        fi
        singularity exec --bind /usr/lib/locale/ {input.img} /opt/deepvariant/bin/run_deepvariant\
            --model_type {params.model} \
            --ref {input.ref} \
            --reads {input.bam} \
            --output_vcf {output} \
            --num_shards {threads} \
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
    input:
        bam = "results/align/raw/{bamtype}/{aligner}_{sample}.aligned.bam",
        img = "/home/jiayiwang/tools/gatk-latest.simg"
    output:
        vcf = "results/variant/GATK/{bamtype}_{aligner}_{sample}.vcf.gz"
    params:
        ref = config["reference_genome"],
        threads = 10,
        pcr_model = choose_pcr_model
    shell:
        """
        tmpbam=$(mktemp --suffix=.bam)
        tmpsrt=$(mktemp --suffix=.bam)
        samtools addreplacerg \
            -r '@RG\\tID:1\\tSM:{wildcards.sample}\\tPL:ONT' \
            -o "$tmpbam" {input.bam}
        
        samtools sort -@ {params.threads} -o "$tmpsrt" "$tmpbam"
        rm "$tmpbam"
        samtools index "$tmpsrt"

        singularity exec --bind /home/jiayiwang/miniconda3/envs/snakemake/bin/python:/usr/bin/python \
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

rule run_longcallR:
    priority: 95
    input: 
        "results/align/raw/{bamtype}/{aligner}_{sample}.aligned.bam"
    output:
        dir=directory("results/variant/longcallR/{bamtype}_{aligner}_{sample}"),
        file="results/variant/longcallR/{bamtype}_{aligner}_{sample}.vcf.gz"
    log:
        "logs/run_longcallR_{bamtype},{aligner},{sample}.log"
    params: 
        threads = 10,
        platform = get_longcallR_platform,
        ref = config["reference_genome"]
    shell:
        """
        if [ ! -f {input}.bai ] || [ {input} -nt {input}.bai ]; then
            samtools index {input};
        fi

        mkdir -p {output.dir}

        ~/tools/longcallR/longcallR \
            --bam-path {input} \
            --ref-path {params.ref} \
            --output {output.dir}/output \
            --preset {params.platform} \
            -t {params.threads} \
            --no-bam-output \
            --min-depth 0 \
            --max-depth 1000000 \
            > {log} 2>&1

        bgzip {output.dir}/output.vcf
        mv {output.dir}/output.vcf.gz {output.file}
        """


CHR = [f"chr{i}" for i in range(1, 23)]

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
    shell:
        """
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
    input:
        vcfs=expand("results/variant/longcallR-nn/{{bamtype}}_{{aligner}}_{{sample}}/vcf/{chr}.vcf", chr=CHR)
    output:
        vcfgz="results/variant/longcallR-nn/{bamtype}_{aligner}_{sample}.vcf.gz",
        tbi="results/variant/longcallR-nn/{bamtype}_{aligner}_{sample}.vcf.gz.tbi"
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


