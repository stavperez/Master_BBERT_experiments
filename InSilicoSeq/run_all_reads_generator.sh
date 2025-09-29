#!/bin/bash
#SBATCH --output=run_all_%j.out
#SBATCH --time=72:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16

# Number of splits to run in parallel
NUM_SPLITS=16

for (( i=0; i<NUM_SPLITS; i++ )); do
    echo "Launching job for split index $i of $NUM_SPLITS"
    # Submit a SLURM job using the run_single.sh script and pass the split parameters.
    sbatch single_run_reads_generator.sh --split-index "$i" --num-splits "$NUM_SPLITS" &
done

wait
echo "All jobs completed."
