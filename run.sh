#!/bin/bash
#SBATCH -p priopark
#SBATCH -A park_contrib
#SBATCH -t 120:00:00
#SBATCH --mem=32G

word=$1
shift

# Always require slurm-drmaa. The SCAN2 commands implicitly require it.
module load slurm-drmaa

# Ensure sentieon is set up
which sentieon
if [ "x$?" == "x1" ]; then
    echo "sentieon is not on your \$PATH and might not be properly enabled. be sure to set the"
    echo "following environment variables:"
    echo "export SENTIEON_LICENSE=address_of_license_server:portnumber"
    echo "export SENTIEON_INSTALL_DIR=/path/to/sentieon_toplevel"
    echo "export PATH=\$PATH:/path/to/sentieon_toplevel/bin"
    exit 1
fi

# $flags is always built added to, never taken from
flags='-s=snakemake/Snakefile --dir=. --latency-wait=60 --resources localjob=1 roadmap_download=10 encode_download=10 ucsc_download=10 --rerun-incomplete'
        #--restart-times 2 \  # This is NECESSARY for some jobs that have step-up memory reqs
jobflag='-j=10'
kgflag=''
drmaaflag=''
usedrmaa='false'

if [ "x$word" == 'xdry' ]; then
    flags="$flags $@ --dryrun --quiet"
elif [ "x$word" == 'xdryreason' ]; then
    flags="$flags $@ --dryrun --reason"
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
    jobflag='-j=6' #20'
    kgflag='--keep-going'
elif [ "x$word" == "xcluster" ]; then
    #echo "be sure to run: module load slurm-drmaa"

    mkdir -p cluster-logs
    usedrmaa='true'
    jobflag='-j=1000'
    kgflag='--keep-going'
    flags="$flags $@ --max-status-checks-per-second=0.5" # --restart-times=2"  # this is handled earlier
    #flags="--max-jobs-per-second 0.05 --max-status-checks-per-second 0.1 --restart-times 2"
else
    echo "usage: $0 {dry,unlock,make_pcawg_metadata,test,local,cluster} [optional snakemake arguments]"
    exit 1
fi


# I can't get $drmaaflags to substitute properly because of the internal 's
if [ $usedrmaa == "true" ]; then
    snakemake $flags $kgflag \
        $jobflag \
        --slurm --default-resources slurm_account=park_contrib slurm_partition=priopark runtime=2880
        #--drmaa=' -p priopark -A park_contrib --mem={resources.mem_mb} -c {threads} -t 4:00:00 -o cluster-logs/slurm-%A.log'
        #--max-threads=12 $jobflag \
else
    echo "snakemake $flags $kgflag --max-threads=12 $jobflag"
    snakemake $flags $kgflag --max-threads=12 $jobflag
fi
