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

cd /projects/intro2gds/I2GDS2025/G2_PlantDisease/Jingjing/results

SPADES_DIR=./Assembled_modern   # modern strain assemble
CHECKM_OUT=./modern_checkm_results   # CheckM output
SUMMARY_FILE=$CHECKM_OUT/checkm_summary.csv
mkdir -p "$CHECKM_OUT"


# ------------------------------
# 3️⃣ loop for every strain 
# ------------------------------
for strain in "$SPADES_DIR"/*; do
    if [ -d "$strain" ]; then
        sample=$(basename "$strain")
        fasta="$strain/scaffolds.fasta"

        # check if fasta exists 
        if [ ! -f "$fasta" ]; then
            echo "No scaffolds.fasta found for $sample, skipping..."
            continue
        fi

        OUTDIR="$CHECKM_OUT/$sample"
        mkdir -p "$OUTDIR"

        echo "Running CheckM for sample: $sample ..."
        checkm lineage_wf \
            -x fasta \
	    --reduced_tree \
            "$strain" \
            "$OUTDIR" \
            --threads 16


    fi
done


# ------------------------------
# 4. Summarize results
# ------------------------------
echo "Generating summary file ..."
checkm qa -o 2 -f "$SUMMARY_FILE" --tab_table "$CHECKM_OUT"/*/lineage.ms "$CHECKM_OUT"

echo "✅ All done!"
echo "Summary saved to: $SUMMARY_FILE"
echo "Job finished at $(date)"


echo "==== All CheckM analyses completed at $(date) ===="
