#!/bin/bash
#SBATCH --job-name=checkm_batch
#SBATCH --account=introtogds
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100GB
#SBATCH --time=48:00:00
#SBATCH --mail-user=jingjingy@vt.edu
#SBATCH --mail-type=ALL

set -o pipefail

echo "==== CheckM batch job started at $(date) ===="


# ------------------------------
# 1️⃣ activate CheckM environment
# ------------------------------
# activate Miniconda
module load Miniconda3/24.7.1-0
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate checkm_env


# ------------------------------
# 2. set work directory and path
# ------------------------------

cd /projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing

# 3. run CheckM
checkm lineage_wf \
    -x fasta \
    --reduced_tree \
    ./spades_output \
    ./results/Historical_checkm_results \
    --threads 16

# generate summary
checkm qa \
    ./results/Historical_checkm_results/lineage.ms \
    ./results/Historical_checkm_results \
    -o 2 \
    -f ./results/Historical_checkm_results/checkm_summary.csv


echo "✅ All done!"
echo "Summary saved to: $SUMMARY_FILE"
echo "Job finished at $(date)"
