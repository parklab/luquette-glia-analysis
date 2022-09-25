#!/bin/bash
#SBATCH -p priopark
#SBATCH -A park_contrib
#SBATCH -t 120:00:00
#SBATCH --mem=32G

word=$1
shift

# $flags is always built added to, never taken from
flags='-s=snakemake/Snakefile --dir=. --latency-wait=60 --rerun-incomplete --rerun-triggers=mtime'
        #--restart-times 2 \  # This is NECESSARY for some jobs that have step-up memory reqs
jobflag='-j=10'
kgflag=''
drmaaflag=''
usedrmaa='false'

if [ "x$word" == 'xdry' ]; then
    flags="$flags $@ --dryrun --quiet" # --reason"
    #flags="$flags --dryrun --reason"
elif [ "x$word" == 'xunlock' ]; then
    flags="$flags --unlock"
elif [ "x$word" == 'xmake_pcawg_metadata' ]; then
    flags="$flags --config make_pcawg_metadata=1 --until metadata/pcawg_metadata.csv"
elif [ "x$word" == 'xtest' ]; then
    flags="$flags $@"
    jobflag='-j=1'
    kgflag=''
elif [ "x$word" == 'xlocal' ]; then
    flags="$flags $@"
    jobflag='-j=8' #20'
    kgflag='--keep-going'
elif [ "x$word" == "xcluster" ]; then
    echo "be sure to run: module load slurm-drmaa"
    usedrmaa='true'
    jobflag='-j=1000'
    kgflag='--keep-going'
    flags="$flags $@ --max-status-checks-per-second=0.1 --restart-times=2"
    #flags="--max-jobs-per-second 0.05 --max-status-checks-per-second 0.1 --restart-times 2"
else
    echo "usage: $0 {dry,unlock,make_pcawg_metadata,test,local,cluster} [optional snakemake arguments]"
    exit 1
fi


# I can't get $drmaaflags to substitute properly because of the internal 's
if [ $usedrmaa == "true" ]; then
    snakemake $flags $kgflag \
        --max-threads=12 $jobflag \
        --drmaa=' -p priopark -A park_contrib --mem={resources.mem} -c {threads} -t 24:00:00 -o cluster-logs/slurm-%A.log'
else
    snakemake $flags $kgflag --max-threads=12 $jobflag
fi
