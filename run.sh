#!/bin/bash
#SBATCH -p priopark
#SBATCH -A park_contrib
#SBATCH -t 120:00:00
#SBATCH --mem=7G

word=$1

jobflag='-j 20'
kgflag=''
flags=''
drmaaflag=''

if [ "x$word" == 'xdry' ]; then
    flags="--dryrun --quiet" # --reason"
elif [ "x$word" == 'xtest' ]; then
    jobflag='-j 1'
    kgflag=''
elif [ "x$word" == 'xlocal' ]; then
    jobflag='-j 20'
    kgflag='--keep-going'
else
    jobflag='-j 1000'
    kgflag='--keep-going'
    drmaaflag="--drmaa ' -p priopark -A park_contrib --mem={resources.mem} -c {threads} -t 24:00:00 -o cluster-logs/slurm-%A.log'"
    flags="--max-status-checks-per-second 0.1" # --restart-times 2"
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
    --max-threads 12 $jobflag $drmaaflag \
    #--max-inventory-time 0 \
