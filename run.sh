#!/bin/bash
#SBATCH -p priopark
#SBATCH -A park_contrib
#SBATCH -t 120:00:00
#SBATCH --mem=7G

dryrun=$1

if [ "x$dryrun" != 'x' ]; then
    flags="--dryrun --quiet" # --reason"
fi

#module load slurm-drmaa

#export TMPDIR=$(realpath try3/tmp)
#echo "TMPDIR=$TMPDIR"

# tmpdir is necessary for using workflow.source to properly find
# scripts from within module calls.
snakemake $flags \
    --dir . \
    --latency-wait 60 \
    --keep-going \
    -s snakemake/Snakefile \
    --restart-times 2 \
    --max-inventory-time 0 \
    --max-threads 12 \
    -j 1000 \
    --drmaa ' -p priopark -A park_contrib --mem={resources.mem} -t 24:00:00 -o cluster-logs/slurm-%A.log'
    #--max-status-checks-per-second 0.1 \
