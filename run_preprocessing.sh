#!/bin/bash

#SBATCH --time=2-12:00:00
#SBATCH --partition=spott
#SBATCH --account=pi-spott
#SBATCH --job-name=pb_preprocessing
#SBATCH --output=%x_%j.out 
#SBATCH --error=%x_%j.err 
#SBATCH --mem=2GB
#SBATCH --ntasks-per-node=1

## Output log and error files will be saved with prefix of the job name and job ID.


###############
##  Parameters : Config file. 
###############
configfile=config_run.yaml

###############
##  MAIN
###############
## Start Run: 
echo Starting Time is `date "+%Y-%m-%d %H:%M:%S"`
start=$(date +%s)


## Run snakemake V7
## Snakemake full path. 
smk=/project/spott/dveracruz/pacbio_preprocessing/workflow/pb_complete.smk

source activate /project/spott/dveracruz/bin/miniconda3/envs/fiber_sq

snakemake -s $smk --configfile $configfile --unlock ## In case it's locked. 

## Produce DAG first. 
snakemake -s $smk --configfile $configfile --dag | dot -Tpng > DAG_run.png

## Main run. 
snakemake -s $smk --configfile $configfile \
    --cluster "sbatch --account=pi-spott --partition={resources.partition} --ntasks-per-node={threads} --mem={resources.mem_mb}" \
    --jobs 50 --keep-target-files --keep-going --rerun-incomplete 

echo Pre-processing: Ending Time is `date "+%Y-%m-%d %H:%M:%S"`
end=$(date +%s)
time=$(( ($end - $start) / 60 ))
echo Used Time is $time mins
