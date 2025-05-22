## Snakemake: PacBio pre-processing pipeline
## Author: Diana Vera Cruz
## Northwestern output data, (HiFi BAM file with modifications called) -> Aligned, no kinetics BAM file.

## This script will run the following steps: 
## 1. Omit Kinetics from raw BAM file.
## 2. Calculate nucleosomes poistions with ft nucleosome.
## 3. pbmm2 align: Align reads to reference genome.

## Running pipeline commands: 
##   tmux new -s this
##   conda activate fiber_sq
## snakemake -s ./workflow/pb_preprocessing.smk  --cluster "sbatch --account=pi-spott --partition=spott --nodes=1 --mem-per-cpu=50000" --jobs 50 --keep-target-files --keep-going --rerun-incomplete
## Dry run: snakemake -s ./workflow/pb_preprocessing.smk  -np --cluster "sbatch --account=pi-spott --partition=caslake --nodes=1 --mem-per-cpu=20000"
## DAG:  snakemake -s ./workflow/pb_preprocessing_v2.smk  --dag | dot -Tpdf > dag_v2.pdf


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
        expand("{dir}/{sample}.5mC.6mA.nuc.aligned.bam", sample = samples.keys(), dir = out_dir),
        expand("{dir}/{sample}.5mC.6mA.nuc.aligned.bam.bai", sample = samples.keys(), dir = out_dir),


## Steps 1&2. 
rule nuc_and_omit_kinetics: 
    input:
        bam = lambda wildcards: samples[wildcards.sample],
    output:
        bam = "{dir}/{sample}.5mC.6mA.nuc.bam"
    threads: 8
    resources:
        mem_mb=80000,
    shell:
        """
        source activate /project/spott/software/conda/pacbio
        ## Omit kinetics with perl script. 
        ft add-nucleosomes {input.bam} - | ft clear-kinetics - {output.bam}
        """ 

## Step 3: Align reads
## Memory and resources: Recommended more memory per thread, than increase threads.
## Updated to 24 threads and 100GB of memory, 10GB sorting, from 12 threads and 80GB of memory + 40GB sorting.
rule align: 
    input:
        bam = "{dir}/{sample}.5mC.6mA.nuc.bam"
    output:
        bam = "{dir}/{sample}.5mC.6mA.nuc.aligned.bam"
    params:
        mmi_ref = "/project/spott/reference/human/GRCh38.p14/hg38p14.mmi"
    threads: 44
    resources:
        mem_mb=240000,
    shell:
        """
        source activate /project/spott/software/conda/pacbio

        ## Set TMPDIR to a scratch directory.
        export TMPDIR=/scratch/midway3/$USER/tmp_pbmm2
        mkdir -p $TMPDIR
        echo $TMPDIR
        
        ## Run alignment:  J 4-> Recommended for Human genome. 
        pbmm2 align {params.mmi_ref} {input.bam} {output.bam} --preset HIFI --sort -j 0 -J 4 --sort-memory 40G --log-level DEBUG
        
        """

## Step 4: Index aligned BAM file
rule index:
    input:
        bam = "{dir}/{sample}.5mC.6mA.nuc.aligned.bam"
    output:
        bam = "{dir}/{sample}.5mC.6mA.nuc.aligned.bam.bai"
    threads: 1
    resources:
        mem_mb = 24000,
    shell:
        """
        module load samtools

        samtools index {input.bam}
        """
