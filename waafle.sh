#!/bin/bash
#SBATCH -c 2                # Number of cores (-c)
#SBATCH -t 3-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --mem=2G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o waafle.out  # File to which STDOUT will be written
#SBATCH -e waafle.err  # File to which STDERR will be written
#SBATCH --mail-type=BEGIN,END,FAIL --mail-user=dichen@coh.org

snakemake -j 50 \
    --cluster "mkdir -p logs/{rule} && sbatch \
        --partition=compute \
        --mem={resources.mem_mb}M \
        --cpus-per-task={threads} \
        --time={resources.runtime} \
        --output=logs/{rule}/{rule}.%j.out \
        --error=logs/{rule}/{rule}.%j.err" \
    --latency-wait 60
