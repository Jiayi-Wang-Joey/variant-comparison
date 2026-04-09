rule preprocess: 
    priority: 100
    input: 
        bam="data/raw/{sample}.bam",
        primer="data/primer/primer_isoseq.fasta"
    output: "data/preprocessed/{sample}.bam"
    conda: "../envs/preprocessing.yaml"
    log:
        stderr="logs/preprocess/{sample}.stderr"
    shell: 
        """
        lima {input.bam} {input.primer} data/raw/{wildcards.sample}.fl.bam --isoseq --peek-guess 2> {log}    
        isoseq refine data/raw/{wildcards.sample}.fl.NEB_5p--NEB_Clontech_3p.bam {input.primer} {output} 2> {log}        
        """    


rule bam2fq:
    input:
        bam="data/{status}/{sample}.bam"
    output:
        "results/align/{status}/bam2fq/{sample}.fastq.gz"
    log:
        stderr="logs/align/{status}/bam2fq/{sample}.stderr"
    threads: config["convert_threads"]
    conda: "../envs/samtools.yaml"
    params:
        threads=config["convert_threads"]
    shell:
        """
        samtools bam2fq -@ {params.threads} {input.bam} | gzip > {output}
        """


def minimap2_preset(wildcards):
    if "IsoSeq" in wildcards.sample:
        return "-ax splice:hq -uf"
    elif "MasSeq" in wildcards.sample:
        return "-ax splice:hq -uf"
    elif "dRNA" in wildcards.sample:
        return "-ax splice -uf -k14"
    elif "cDNA" in wildcards.sample:
        return "-ax splice"
    else:
        raise ValueError(f"Cannot determine minimap2 preset from sample name: {wildcards.sample}")


# rule minimap2_align:
#     priority: 98
#     input:
#         reads = rules.bam2fq.output,
#         transcriptome = config["transcriptome_bed"],
#         genome = config["reference_genome"],
#     output:
#         "results/align/{status}/origin/minimap2_{sample}.aligned.bam"
#     params:
#         align_map_bam_threads = config["align_map_bam_threads"],
#         align_sort_bam_threads = config["align_sort_bam_threads"],
#         align_sort_bam_memory_gb = config["align_sort_bam_memory_gb"],
#         preset = minimap2_preset
#     conda:
#         "../envs/minimap2.yaml"
#     log:
#         stdout = "logs/minimap2/{status}/{sample}.out",
#         stderr = "logs/minimap2/{status}/{sample}.err"
#     # wildcard_constraints:
#     #     status="raw|preprocessed" if "{sample}" in ISOSEQ_SAMPLE else "raw" 
#     shell:
#         """
#         minimap2 {params.preset} --junc-bed {input.transcriptome} \
#             -t {params.align_map_bam_threads} \
#             {input.genome} {input.reads} | samtools sort \
#             -@ {params.align_sort_bam_threads} \
#             -m{params.align_sort_bam_memory_gb}g \
#             -o {output} > {log.stdout} 2> {log.stderr}
#         """


# rule pbmm2_align:
#     priority: 98
#     input:
#         reads = "results/align/{status}/bam2fq/{sample}.fastq.gz",   
#         genome = config["reference_genome"]
#     output:
#         bam = "results/align/{status}/origin/pbmm2_{sample}.aligned.bam"
#     params:
#         preset = "ISOSEQ",  
#         threads = config["align_map_bam_threads"]
#     conda:
#         "../envs/pbmm2.yaml"
#     log:
#         stdout = "logs/pbmm2/{status}/{sample}.out",
#         stderr = "logs/pbmm2/{status}/{sample}.err"
#     # wildcard_constraints:
#     #     status="raw|preprocessed" if "{sample}" in ISOSEQ_SAMPLE else "raw" 
#     shell:
#         """
#         pbmm2 align \
#             --preset {params.preset} \
#             --sort \
#             --num-threads {params.threads} \
#             {input.genome} \
#             {input.reads} \
#             {output.bam} > {log.stdout} 2> {log.stderr}
#         """

# rule split_ncigar:
#     priority: 97
#     input:
#         bam = "results/align/{status}/origin/{aligner}_{sample}.aligned.bam",
#         ref = config["reference_genome"],
#         img = "/home/jiayiwang/tools/gatk-latest.simg"  
#     output: 
#         "results/align/{status}/splitNC/{aligner}_{sample}.aligned.bam"
#     log:
#         "logs/split_ncigar_reads_{status}_{aligner}_{sample}.log"
#     params:
#         threads = 10
#     wildcard_constraints:
#         status="raw|preprocessed" if "{sample}" in ISOSEQ_SAMPLE else "raw" 
#     shell:
#         """
#         singularity exec --bind /home/jiayiwang/miniconda3/envs/snakemake/bin/python:/usr/bin/python {input.img} /gatk/gatk --java-options \
#         "-Xmx16G -XX:+UseParallelGC -XX:ParallelGCThreads={params.threads} -Djava.util.concurrent.ForkJoinPool.common.parallelism={params.threads} -Djava.io.tmpdir=/home/jiayiwang/tmp1" SplitNCigarReads \
#             -R {input.ref} \
#             -I {input.bam} \
#             -O {output} \
#         > {log} 2>&1
#         """ 

# rule flag_correction:
#     priority: 96
#     input:
#         bam = "results/align/{status}/origin/{aligner}_{sample}.aligned.bam",
#         sncr_bam = rules.split_ncigar.output
#     output: 
#         "results/align/{status}/transformed/{aligner}_{sample}.aligned.bam"
#     log:
#         "logs/flag_correction_{status}_{aligner}_{sample}.log"
#     params:
#         threads = 5
#     wildcard_constraints:
#         status="raw|preprocessed" if "{sample}" in ISOSEQ_SAMPLE else "raw" 
#     shell:
#         """
#         Rscript /home/jiayiwang/tools/lrRNAseqVariantCalling/tools/flagCorrection.r \
#             {input.bam} \
#             {input.sncr_bam} \
#             {output} \
#             {params.threads} \
#         > {log} 2>&1
#         rm {input.sncr_bam}
#         """


rule bam_quality:
    priority: 80
    input:
        "results/align/{status}/origin/{aligner}_{sample}.aligned.bam"
    output:
        stat = "results/align/{status}/quality/{aligner}_{sample}.tsv",
        cov = "results/align/{status}/coverage/{aligner}_{sample}.tsv"
    log:
        stdout = "logs/quality/{status}_{aligner}_{sample}.log",
        stderr = "logs/quality/{status}_{aligner}_{sample}.err"
    shell:
        """
        echo -e "sample\\tstatus\\taligner\\tmetric\\tvalue\\tannotation" > {output.stat}
        echo -e "sample\\tstatus\\taligner\\trange\\tdepth\\tcount" > {output.cov}

        samtools stats {input} 2>> {log.stderr} | tee >(grep ^SN | \
            awk -v sample="{wildcards.sample}" -v status="{wildcards.status}" -v aligner="{wildcards.aligner}" 'BEGIN {{FS="\\t"; OFS="\\t"}} \
            {{metric = $2; value = $3; annotation = ($0 ~ /#/) ? substr($0, index($0, "#")) : ""; print sample, status, aligner, metric, value, annotation}}' >> {output.stat}) \
        | grep ^COV | cut -f 2- | \
            awk -v s="{wildcards.sample}" -v t="{wildcards.status}" -v a="{wildcards.aligner}" 'BEGIN {{OFS="\\t"}} {{print s,t,a,$1,$2,$3}}' >> {output.cov}
        """