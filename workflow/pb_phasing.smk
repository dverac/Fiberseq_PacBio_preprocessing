## Snakemake: PacBio Phasing pipeline
## Author: Diana Vera Cruz
## Uses Aligned BAM files from PacBio and performs phasing using HiPhase.

## This script will run the following steps: 
## 1. Subset reference VCF: Subset reference VCF to sample name. (BCFtools)
## 2. Rename SM tag: Rename SM tag in BAM file to match sample name in VCF. (Samtools)
## 3. Phasing: Phasing of reads using HiPhase. (HiPhase)

## Running pipeline commands:
##   tmux new -s this
##   conda activate fiber_sq
## snakemake -s ./workflow/pb_phasing.smk  --cluster "sbatch --account=pi-spott --partition=spott --nodes=1 --mem-per-cpu=50000" --jobs 50 --keep-target-files --keep-going --rerun-incomplete
## Dry run: snakemake -s ./workflow/pb_phasing.smk  -np --cluster "sbatch --account=pi-spott --partition=caslake --nodes=1 --mem-per-cpu=20000"
## DAG:  snakemake -s ./workflow/pb_phasing.smk  --dag | dot -Tpdf > dag_phasing.pdf


##########
## Samples PARAMETERS
## Read config file: Input and output directories. 
# configfile: "config_run_complete.yaml" -> Will be produced in the sbatch file that runs the pipeline.(workflow/generate_config_run.r)
##########
configfile: "/project/spott/dveracruz/pacbio_preprocessing/workflow/config_run_complete.yaml"  

samples = config['samples'] ## Dictionary with sample names and paths to raw bam files.
out_dir = config['out_dir'] ## Output directory
vcf = config['vcf'] ## Reference VCF file for phasing.
vcf_sm = config['vcf_sm'] ## Sample name (SM tag) in VCF file per sample.


##########
## RULES
###########

## Rule all: 
rule all: 
    input:
        expand("{dir}/{sample}.5mC.6mA.aligned.phased.bam", sample = samples.keys(), dir = out_dir),

## Step 1: Subset reference VCF to sample name.
rule subset_vcf:
    input:
        bam = "{dir}/{sample}.5mC.6mA.aligned.bam",
        ref_vcf = "/project/spott/reference/YRI_genotypes/LCL.vcf.gz",
    output:
        vcf = "{dir}/{sample}/{sample}.vcf.gz"
    params:
        sample = "{sample}",
        ref_fa = "/project/spott/reference/human/GRCh38.p14/hg38.fa",
        bcftools = "/project/spott/software/bcftools",
    threads: 24
    resources:
        mem_mb=200000,
        partition="spott",
    shell:
        """
        {params.bcftools}/bcftools view -Oz -s {params.sample} {input.ref_vcf} > {output.vcf}
        {params.bcftools}/bcftools index {output.vcf}
        """

## Correct Sample names: Match SM tag to sample name in VCF to be used for phasing.
rule SMcorrect:
    input:
        bam = lambda wildcards: samples[wildcards.sample]
    output:
        bam = "{dir}/{sample}.smcorrect.bam"
    params:
        sm = lambda wildcards: vcf_sm[wildcards.sample],
    threads: 2
    resources:
        mem_mb=10000,
    shell:
        """
        module load samtools
        mkdir -p $(dirname {output.bam})

        samtools view -H {input.bam}  | sed "s/SM:[^\t]*/SM:{params.sm}/g" | samtools reheader - {input.bam} > {output.bam}
        samtools index {output.bam}
        """
        

## Functional part, ref VCF is LC_18508.vcf.gz
## Low coverage: 4 threads, 20GB
## High coverage: 8 threads, 60GB
rule phasing:
    input:
        bam = "{dir}/{sample}.smcorrect.bam",
        vcf = config["vcf"], 
    output:
        bam = "{dir}/{sample}.5mC.6mA.aligned.phased.bam",
        vcf = "{dir}/{sample}/{sample}.5mC.6mA.aligned.phased.vcf.gz",
        haplotag = "{dir}/{sample}/read-level-phasing.tsv",
        summary = "{dir}/{sample}/summary.tsv",
        stats = "{dir}/{sample}/stats.tsv",
        blocks = "{dir}/{sample}/blocks.tsv",
    params:
        ref_fa = "/project/spott/reference/human/GRCh38.p14/hg38.fa",
        hiphase = "/project/spott/software/hiphase_v1.4.5/hiphase-v1.4.5-x86_64-unknown-linux-gnu",
        vcf_sm = lambda wildcards: vcf_sm[wildcards.sample],
    threads: 8
    resources:
        mem_mb=60000,
    shell:
        """
        {params.hiphase}/hiphase --threads {threads} \
            --ignore-read-groups \
            --reference {params.ref_fa} \
            --bam {input.bam} \
            --vcf {input.vcf} \
            --sample-name {params.vcf_sm} \
            --output-bam {output.bam} \
            --output-vcf {output.vcf} \
            --haplotag-file {output.haplotag} \
            --summary-file {output.summary} \
            --stats-file {output.stats} \
            --blocks-file {output.blocks}
        """



