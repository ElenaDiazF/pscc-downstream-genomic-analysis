#!/bin/bash
#SBATCH --job-name=pscc_reports
#SBATCH --partition=long
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=9
#SBATCH --mem=150G
#SBATCH --output=output/logs/reports_%j.out
#SBATCH --error=output/logs/reports_%j.err

# Print job information
echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "Start Time: $(date)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: ${SLURM_MEM_PER_NODE}MB"
echo "================================================="
echo ""

# ============================================================================
# Setup
# ============================================================================
source /PROJECTES/SQUAMOLAB/pipelines/WES2/scripts/lib.sh

# Load R environment
echo "⏳ Loading R environment (r-env)..."
conda_activate r-env
echo ""

# Strict error handling (after conda activation to avoid unbound variable errors in conda scripts)
set -euo pipefail

# Navigate to the working directory
cd /PROJECTES/SQUAMOLAB/projects/pscc/wes
# Create log directory
mkdir -p output/logs

# Run the R script
echo "Starting R script execution..."
Rscript output/reports_list12.R

# Check exit status
EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "=================================================="
    echo "Job completed successfully!"
    echo "End Time: $(date)"
    echo "=================================================="
else
    echo ""
    echo "=================================================="
    echo "Job failed with exit code: $EXIT_CODE"
    echo "End Time: $(date)"
    echo "=================================================="
    exit $EXIT_CODE
fi
