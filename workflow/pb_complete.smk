## Snakemake: PacBio pre-processing pipeline
## Author: Diana Vera Cruz
## Latest update: 04/17/2025
## Requires data with 5mC and 6mA methylation information: Output from Northwestern since early 2025. 
## Process (HiFi BAM file) with MM/ML tags to create aligned, phased BAM file with methylation information for 5mC and 6mA.
## Requires reference genome and VCF file for phasing per sample. 

## This script will run the following steps: 
## 1&2. ft add-nucleosomes + ft clear-kinetics: Add nucleosome information per read, and clear kinetics to reduce size. 
## 3. pbmm2 align: Align reads to reference genome, sort and index BAM file.
## 4&5. HiPhase: SM correctiong and Phasing (Optional)

## Running pipeline commands: 
##   tmux new -s this
##   conda activate fiber_sq
## snakemake -s ./workflow/pb_complete.smk --configfile ./workflow/config_run_complete.yaml  --cluster "sbatch --account=pi-spott --partition={resources.partition} --ntasks-per-node={threads} --mem={resources.mem_mb}"--jobs 50 --keep-target-files --keep-going --rerun-incomplete
## Dry run: snakemake -s ./workflow/pb_complete.smk --configfile ./workflow/config_run_complete.yaml  -np 
## DAG:  snakemake -s ./workflow/pb_complete.smk --configfile ./workflow/config_run_complete.yaml --dag | dot -Tpdf > dag.pdf


##########
## Samples PARAMETERS
## Read config file: Input and output directories. 
## Default config file: "/project/spott/dveracruz/pacbio_preprocessing/config_run.yaml"
# configfile: "config_run.yaml" -> Will be produced in the sbatch file that runs the pipeline.(workflow/generate_config_run.r)
##########

## 1) First check that the elements in loaded config are correct. 

## Check config contains samples, out_dir and phasing.
for key in ['samples', 'out_dir', 'phasing']:
    value = config.get(key)
    if not value:
        raise ValueError(f"Missing or empty required config key: '{key}'")

samples = config['samples'] ## Dictionary with sample names and paths to raw bam files.
out_dir = config['out_dir'] ## Output directory
phasing = config['phasing'] ## Should we phase the reads?

## 2) Check output directory exists, if not create it. 
import os
os.makedirs(out_dir, exist_ok=True)

## 3) If phasing is yes, check that VCF and vcf_sm are present, otherwise, print, no phasing. 
if phasing.lower() == "yes":
    ## Check vcf and vcf_sm elements exist.
    for key in ['vcf', 'vcf_sm']:
        value = config.get(key)
        if not value:
            raise ValueError(f"Missing or empty required config key when phasing = 'Yes': '{key}'")

    ## Check the VCF exists. 
    if not os.path.exists(config['vcf']):
        raise FileNotFoundError(f"Reference VCF file does not exist: {config['vcf']}")

    ## Check that all the keys in samples dictionary are present in the vcf_sm dictionary.
    vcf_sm = config['vcf_sm']
    missing_keys = set(samples.keys()) - set(vcf_sm.keys())
    if missing_keys:
        raise KeyError(f"The following keys in the samples dictionary are missing in the vcf_sm dictionary: {missing_keys}")
else:
    print("Phasing is set to No. Skipping phasing step.", file=sys.stderr)

##############
## Resources style: Set up different resources for high-coverage samples versus low coverage, specially for alignment & phasing. 
if 'high_cov' not in config:
    config['high_cov'] = 'no' ## Default to no high coverage.
    
if(config['high_cov'].lower() == 'yes'):
    print("High coverage samples, setting up resources for high coverage.", file=sys.stderr)
    ## Make a dictionary with threads and memory for each sample.
    ## align: pbmm2 --sort -j 40 -J 4 --sort-memory 40G
    res = {
        "align_partition": "bigmem", ## Use bigmem. 
        "align_threads": 44, 
        "align_mem_mb": 240000,
        "align_sort_mem": 40
    }
else:
    print("Low coverage samples, setting up resources for low coverage.", file=sys.stderr)
    ## align: pbmm2 --sort -j 28 -J 4 --sort-memory 10G
    res = {
        "align_partition": "caslake", ## Use caslake.
        "align_threads":32, 
        "align_mem_mb": 120000, ## Updated to 120GB of memory to avoid out of memory errors.
        "align_sort_mem": 10
    }


##########
## RULES
###########

## Rule all: If phasing is Yes, add phased bam file to output. 
phasing_output = (
    expand("{dir}/{sample}.5mC.6mA.aligned.phased.bam", sample=samples.keys(), dir=out_dir)
    if phasing.lower() == "yes"
    else []
)

rule all: 
    input:
        expand("{dir}/{sample}.5mC.6mA.nuc.aligned.bam.bai", sample = samples.keys(), dir = out_dir),
        *phasing_output


## Steps 1&2. 
rule nuc_and_omit_kinetics: 
    input:
        bam = lambda wildcards: samples[wildcards.sample],
    output:
        bam = "{dir}/{sample}.5mC.6mA.nuc.bam"
    threads: 8
    resources:
        mem_mb=80000,
        partition="caslake",
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
        mmi_ref = "/project/spott/reference/human/GRCh38.p14/hg38p14.mmi",
        sort_mem = res["align_sort_mem"],
    threads: res["align_threads"]
    resources:
        mem_mb = res["align_mem_mb"],
        partition = res["align_partition"],
    shell:
        """
        source activate /project/spott/software/conda/pacbio

        ## Set TMPDIR to a scratch directory.
        export TMPDIR=/scratch/midway3/$USER/tmp_pbmm2
        mkdir -p $TMPDIR
        echo $TMPDIR
        
        ## Run alignment:  J 4-> Recommended for Human genome. 4 threads for sorting.
        pbmm2 align {params.mmi_ref} {input.bam} {output.bam} --preset HIFI --sort -j 0 -J 4 --sort-memory {params.sort_mem}G --log-level DEBUG
        
        """

## Step 3B: Index aligned BAM file
rule index:
    input:
        bam = "{dir}/{sample}.5mC.6mA.nuc.aligned.bam"
    output:
        bam = "{dir}/{sample}.5mC.6mA.nuc.aligned.bam.bai"
    threads: 1
    resources:
        mem_mb = 32000,
        partition = 'caslake',
    shell:
        """
        module load samtools

        samtools index {input.bam}
        """


## Step 4 & 5: CorrectSM & Phasing (Optional)
if phasing == "Yes":
    ## Step 4: Correct SM Sample names: Match SM tag to sample name in VCF to be used for phasing.
    rule SMcorrect:
        input:
            bam = "{dir}/{sample}.5mC.6mA.nuc.aligned.bam"
        output:
            bam = "{dir}/{sample}.smcorrect.bam"
        params:
            sm = lambda wildcards: vcf_sm[wildcards.sample],
        threads: 2
        resources:
            mem_mb=10000,
            partition="caslake",
        shell:
            """
            module load samtools
            mkdir -p $(dirname {output.bam})

            samtools view -H {input.bam}  | sed "s/SM:[^\t]*/SM:{params.sm}/g" | samtools reheader - {input.bam} > {output.bam}
            samtools index {output.bam}
            """
    ## Step 5: Phasing
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
            partition="caslake",
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

   

