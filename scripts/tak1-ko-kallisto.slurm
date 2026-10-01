#!/bin/bash

#SBATCH --job-name=tak1-ko-kallisto
#SBATCH --output=/work/newton_lab/ag_analysis/a549-tak1-ko/logs/%j_kallisto.out
#SBATCH --error=/work/newton_lab/ag_analysis/a549-tak1-ko/logs/%j_kallisto.err
#SBATCH --time=02:00:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1

module load kallisto/0.46.1
cd /work/newton_lab/ag_analysis/a549-tak1-ko/

index=/work/newton_lab/ag_analysis/ref_seq/hg38_p14_kallisto.idx
fastq1=$1
fastq2=$2
output=$3
mkdir -p "$(dirname "$output")"

kallisto quant \
    -i "$index" \
    -o "$output" \
    --threads 4 \
    -b 100 \
    "$fastq1" "$fastq2"
