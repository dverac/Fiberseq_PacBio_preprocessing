## Snakemake: PacBio pre-processing pipeline
## Author: Diana Vera Cruz
## Process raw data (HiFi BAM file) to create methylated calls and aligned BAM file. 

## This script will run the following steps: 
## 1. Jasmine: Call 5mC methylation in raw bam. 
## 2. ft predict-m6a: Call 6mA methylation in 5mC bam. 
## 3. pbmm2 align: Align reads to reference genome.

## Running pipeline commands: 
##   tmux new -s this
##   conda activate fiber_sq
## snakemake -s ./workflow/pb_preprocessing.smk  --cluster "sbatch --account=pi-spott --partition=spott --nodes=1 --mem-per-cpu=50000" --jobs 50 --keep-target-files --keep-going --rerun-incomplete
## Dry run: snakemake -s ./workflow/pb_preprocessing.smk  -np --cluster "sbatch --account=pi-spott --partition=caslake --nodes=1 --mem-per-cpu=20000"
## DAG:  snakemake -s ./workflow/pb_preprocessing.smk  --dag | dot -Tpdf > dag.pdf


##########
## Samples PARAMETERS
## Read config file: Input and output directories. 
# configfile: "config_run.yaml" -> Will be produced in the sbatch file that runs the pipeline.(workflow/generate_config_run.r)
##########
configfile: "config_run.yaml"  ## Default, 

samples = config['samples'] ## Dictionary with sample names and paths to raw bam files.
out_dir = config['out_dir'] ## Output directory

##########
## RULES
###########

## Rule all: 
rule all: 
    input:
        expand("{dir}/{sample}.5mC.6mA.aligned.bam", sample = samples.keys(), dir = out_dir)

## Step 1. Jasmine: Call 5mC methylation. 
rule jasmine_5mC: 
    input:
        bam = lambda wildcards: samples[wildcards.sample],
    output:
        bam = "{dir}/{sample}.5mC.bam"
    threads: 12
    resources:
        mem_mb=50000,
    shell:
        """
        module load gcc/13.2.0
        module load cmake/3.26
        module load python/anaconda-2022.05
        source activate /project/spott/software/conda/pacbio

        ## Run Jasmine: Keep kinetics for 6mA calling.
        jasmine {input.bam} {output.bam} -j 32 --log-level INFO --keep-kinetics
        """ 

## Step 2: 6mA methylation with FiberSeq
rule fiberseq_6mA: 
    input:
        bam = "{dir}/{sample}.5mC.bam"
    output:
        bam = "{dir}/{sample}.5mC.6mA.bam"
    threads: 12
    resources:
        mem_mb=50000,
    shell:
        """
        module load gcc/13.2.0
        module load cmake/3.26
        module load python/anaconda-2022.05
        #source activate /project/spott/dveracruz/bin/miniconda3/envs/fiber_sq
        source activate /project/spott/software/conda/pacbio

        ft predict-m6a -t 32 {input.bam} {output.bam}

        source activate /project/spott/software/conda/pacbio

        pbindex {output.bam}
        """ 

## Step 3: Align reads
rule align: 
    input:
        bam = "{dir}/{sample}.5mC.6mA.bam"
    output:
        bam = "{dir}/{sample}.5mC.6mA.aligned.bam"
    params:
        mmi_ref = "/project/spott/reference/human/GRCh38.p14/hg38p14.mmi"
    threads: 12
    resources:
        mem_mb=50000,
    shell:
        """
        module load gcc/13.2.0
        module load cmake/3.26
        module load python/anaconda-2022.05
        source activate /project/spott/software/conda/pacbio

        pbmm2 align {params.mmi_ref} {input.bam} {output.bam} --preset HIFI --sort -j 0 -J 4 --sort-memory 10G --log-level DEBUG
        """

