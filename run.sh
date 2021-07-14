#!/bin/bash
#SBATCH -p priopark
#SBATCH -A park_contrib
#SBATCH -t 120:00:00
#SBATCH --mem=7G

dryrun=$1

if [ "x$dryrun" != 'x' ]; then
    flags="--dryrun --quiet"
fi

module load slurm-drmaa

#export TMPDIR=$(realpath try3/tmp)
#echo "TMPDIR=$TMPDIR"

# tmpdir is necessary for using workflow.source to properly find
# scripts from within module calls.
snakemake $flags \
    --dir try3 \
    --latency-wait 60 \
    --keep-going \
    -s Snakefile \
    --max-inventory-time 0 \
    -j 1000 \
    --max-status-checks-per-second 0.1 \
    --drmaa ' -p park -A park_contrib --mem={resources.mem} -t 12:00:00 -o /n/data1/hms/dbmi/park/jluquette/glia/analysis/try3/cluster-logs/slurm-%A.log'
