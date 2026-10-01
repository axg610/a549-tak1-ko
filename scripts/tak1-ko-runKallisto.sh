# feeds the fastq1, fastq2, and kallisto_path columns of meta file into the slurm job
# sends one job for each row of the meta file.

awk -F'\t' 'NR>1 {
  cmd = "sbatch tak1-ko-kallisto.slurm " \
        $8 " " \
        $9 " " \
        $10
  system(cmd)
}' tak1-ko-meta.txt