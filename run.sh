#!/bin/bash
#SBATCH -p priopark
#SBATCH -A park_contrib
#SBATCH -t 120:00:00
#SBATCH --mem=7G

word=$1

jobflag='-j 20'
kgflag=''
flags=''
drmaaflag="--drmaa \' -p priopark -A park_contrib --mem={resources.mem} -c {threads} -t 24:00:00 -o cluster-logs/slurm-%A.log\'"

if [ "x$word" == 'dry' ]; then
    flags="--dryrun --quiet" # --reason"
elif [ "x$word" == 'test' ]; then
    jobflag='-j 1'
    kgflag=''
else
    kgflag='--keep-going'
fi

#module load slurm-drmaa

#export TMPDIR=$(realpath try3/tmp)
#echo "TMPDIR=$TMPDIR"

# tmpdir is necessary for using workflow.source to properly find
# scripts from within module calls.
snakemake $flags \
    --dir . \
    --latency-wait 60 $kgflag \
    -s snakemake/Snakefile \
    --max-inventory-time 0 \
    --max-threads 12 $jobflag $drmaaflag \
    #--restart-times 2 \
    #--max-status-checks-per-second 0.1 \
