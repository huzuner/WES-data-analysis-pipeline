#WHOLE EXOME SEQUENCING VARIANT CALLING WORKFLOW - from Mapping to Variant calling#
#by WeepingM#

configfile: "config.yaml"

#The following commented part creates a folder named recal and creates empty .recalibrated.bam files for a given number of patients with 1 tumor and 1 blood samples.

###########################################################

#UNCOMMENT THIS PART IN THE FIRST (1ST) RUN#

#os.mkdir("recal")

#for i in range (1,config["n_of_patients"]):
#    os.system ("touch "+"recal/Patient"+str(i)+"_"+"Tumor1_recalibrated.bam")
#    os.system ("touch "+"recal/Patient"+str(i)+"_"+"Blood1_recalibrated.bam")

#UNCOMMENT THIS PART IN THE FIRST (1ST) RUN#

###########################################################

rule all:
    input:
        expand("mapped/{sample}.sam", sample=config["samples"]),
        expand("sorted/{sample}.bam", sample=config["samples"]),
        expand("mapped/stat/{sample}_stat.txt", sample=config["samples"]),
        expand("fixed-rg/{sample}.bam", sample=config["samples"]),
        expand("dedupped/{sample}.bam", sample=config["samples"]),
        expand("dedupped/{sample}.metrics.txt", sample=config["samples"]),
        expand("recal/{sample}_recal_data.table", sample=config["samples"]),
        expand("recal/{sample}_recalibrated.bam", sample=config["samples"]),
        expand("variant_call/{patient}_Tumor_Normal.raw.vcf.gz", patient=config["patients"])
        
rule bwa_map:
    input:
        ref=config["genome"],
        sample=lambda wildcards: expand(f"{config['samples'][wildcards.sample]}_{{num}}.fastq", num=[1,2])
    output:
        "mapped/{sample}.sam"
    threads:
        config["bwa_mem"]
    log:
        "logs/bwa_mem/{sample}.log"
    shell:
        "bwa mem -t {threads} {input} > {output}"

rule samtools_sort:
    input:
        "mapped/{sample}.sam"
    output:
        "sorted/{sample}.bam"
    log:
        "logs/samtools/sort/{sample}.log"
    shell:
        "samtools sort -@ {threads} -o {output} {input}"

rule align_stat:
    input:
        "sorted/{sample}.bam"
    output:
        "mapped/stat/{sample}_stat.txt"
    log:
        "logs/samtools/flagstat/{sample}.log"
    shell:
        "samtools flagstat {input} > {output}"

#POST-ALIGNMENT PROCESSING

rule add_rg:
    input:
        "sorted/{sample}.bam"
    output:
        "fixed-rg/{sample}.bam"
    log:
        "logs/picard/replace_rg/{sample}.log"
    params:
        general="RGLB={sample}_lib RGPL=illumina RGPU={sample}_unit RGSM={sample} VALIDATION_STRINGENCY=LENIENT",
        picard=config["PathToPicard"]
    shell:
        "java -jar {params.picard}/picard.jar AddOrReplaceReadGroups I={input} O={output} {params.general}"

rule mark_dup:
    input:
        "fixed-rg/{sample}.bam"
    output:
        bam="dedupped/{sample}.bam",
        metrics="dedupped/{sample}.metrics.txt"
    log:
        "logs/picard/mark_dup/{sample}.log"
    params:
        general="CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT",
        picard=config["PathToPicard"]
    shell:
        "java -jar {params.picard}/picard.jar MarkDuplicates I={input} O={output.bam} M={output.metrics} {params.general}"

#BASE RECALIBRATION
rule base_recal: #includes both the base recalibrator and apply bqsr steps
    input:
        bam="dedupped/{sample}.bam",
        ref="data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        known1=config["Mills"],
        known2=config["dbSNP"]
    output:
        bam="recal/{sample}_recal_data.table"
    log:
        "logs/gatk/bqsr/base_recal/{sample}.log"
    shell:
        "gatk BaseRecalibrator -I {input.bam} -R {input.ref} --known-sites {input.known1}"
        " --known-sites {input.known2} -O {output}"

rule appy_bqsr:
    input:
        bam="dedupped/{sample}.bam",
        ref="data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        recal_file="recal/{sample}_recal_data.table",
    output:
        bam="recal/{sample}_recalibrated.bam"
    log:
        "logs/gatk/bqsr/apply_bqsr/{sample}.log"
    shell:
        "gatk ApplyBQSR -R {input.ref} -I {input.bam} --bqsr-recal-file {input.recal_file} -O {output}"


rule variant_call:
    input:
        ref=config["genome"],
        tumor="recal/{patient}_Tumor1_recalibrated.bam",
        normal="recal/{patient}_Blood1_recalibrated.bam"
    output:
        "variant_call/{patient}_Tumor_Normal.raw.vcf.gz"
    params:
        gnomad=config["germline"],
        java=config["mutect2_java_opt"]
    log:
        "logs/gatk/mutect2/{patient}.log"
    threads: config["mutect2"]
    shell:
        " gatk {params.java} Mutect2 -R {input.ref} -I {input.normal} -I {input.tumor} -normal '{wildcards.patient}_Blood1'"
        " -O {output} --independent-mates --germline-resource {params.gnomad}"
