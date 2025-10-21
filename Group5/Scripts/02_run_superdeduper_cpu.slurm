#!/bin/bash
#SBATCH --job-name=sdeduper_cpu
#SBATCH --account=introtogds                             # allocation name 
#SBATCH --partition=normal_q
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -D /projects/intro2gds/I2GDS2025/G5_MG_AMR/after_deduper
#SBATCH --output=/projects/intro2gds/I2GDS2025/G5_MG_AMR/logs/%x_%j.out
#SBATCH --error=/projects/intro2gds/I2GDS2025/G5_MG_AMR/logs/%x_%j.err

set -euxo pipefail
trap 'echo "[FATAL] Script failed at line $LINENO"; env | sort | grep -E "SLURM|CONDA|PATH" || true' ERR

INDIR="/projects/intro2gds/I2GDS2025/G5_MG_AMR/Trim_Galore"
OUTDIR="/projects/intro2gds/I2GDS2025/G5_MG_AMR/after_deduper"
LOGDIR="/projects/intro2gds/I2GDS2025/G5_MG_AMR/logs"
mkdir -p "$OUTDIR" "$LOGDIR"

echo "[INFO] WorkDir: $(pwd)"
echo "[INFO] INDIR:   $INDIR"
echo "[INFO] OUTDIR:  $OUTDIR"
echo "[INFO] LOGDIR:  $LOGDIR"

# Activate  environment                               
source /projects/intro2gds/I2GDS2025/tools/miniconda3/etc/profile.d/conda.sh
conda activate htstream12
command -v hts_SuperDeduper

# --------- Collect R1 List (Covering Trim Galore Naming)--------
#   *_1_val_1.fq(.gz) / *_2_val_2.fq(.gz)       
#   *_R1.fq(.gz)   / *_R2.fq(.gz)
#   *_1.fq(.gz)    / *_2.fq(.gz)
#   *_R1.fastq(.gz)/ *_R2.fastq(.gz)
#   *_1.fastq(.gz) / *_2.fastq(.gz)

shopt -s nullglob

mapfile -t R1_LIST < <(find "$INDIR" -maxdepth 1 -type f \( \
   -name "*_1_val_1.fq.gz"    -o -name "*_1_val_1.fastq.gz"    -o \
   -name "*_1_val_1.fq"       -o -name "*_1_val_1.fastq"       -o \
   -name "*_R1.fastq.gz"      -o -name "*_R1.fq.gz"            -o \
   -name "*_1.fastq.gz"       -o -name "*_1.fq.gz"             -o \
   -name "*_R1.fastq"         -o -name "*_R1.fq"               -o \
   -name "*_1.fastq"          -o -name "*_1.fq" \
   \) | sort)

if [[ ${#R1_LIST[@]} -eq 0 ]]; then
  echo "[ERROR] No R1 files found in $INDIR (please check file extensions and naming)"
  exit 2
fi

pairs_done=0
for R1 in "${R1_LIST[@]}"; do
  base="$(basename "$R1")"
  R2=""
  SAMPLE=""

  case "$base" in
    *_1_val_1.fq.gz)
      R2="${R1/_1_val_1.fq.gz/_2_val_2.fq.gz}"
      SAMPLE="${base%_1_val_1.fq.gz}"
      ;;
    *_1_val_1.fastq.gz)
      R2="${R1/_1_val_1.fastq.gz/_2_val_2.fastq.gz}"
      SAMPLE="${base%_1_val_1.fastq.gz}"
      ;;
    *_1_val_1.fq)
      R2="${R1/_1_val_1.fq/_2_val_2.fq}"
      SAMPLE="${base%_1_val_1.fq}"
      ;;
    *_1_val_1.fastq)
      R2="${R1/_1_val_1.fastq/_2_val_2.fastq}"
      SAMPLE="${base%_1_val_1.fastq}"
      ;;
    *_R1.fastq.gz)
      R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
      SAMPLE="${base%_R1.fastq.gz}"
      ;;
    *_R1.fq.gz)
      R2="${R1/_R1.fq.gz/_R2.fq.gz}"
      SAMPLE="${base%_R1.fq.gz}"
      ;;
    *_1.fastq.gz)
      R2="${R1/_1.fastq.gz/_2.fastq.gz}"
      SAMPLE="${base%_1.fastq.gz}"
      ;;
    *_1.fq.gz)
      R2="${R1/_1.fq.gz/_2.fq.gz}"
      SAMPLE="${base%_1.fq.gz}"
      ;;
    *_R1.fastq)
      R2="${R1/_R1.fastq/_R2.fastq}"
      SAMPLE="${base%_R1.fastq}"
      ;;
    *_R1.fq)
      R2="${R1/_R1.fq/_R2.fq}"
      SAMPLE="${base%_R1.fq}"
      ;;
    *_1.fastq)
      R2="${R1/_1.fastq/_2.fastq}"
      SAMPLE="${base%_1.fastq}"
      ;;
    *_1.fq)
      R2="${R1/_1.fq/_2.fq}"
      SAMPLE="${base%_1.fq}"
      ;;
    *)
      echo "[WARN] Unmatched naming pattern:$base"; continue;;
  esac

  if [[ ! -f "$R2" ]]; then
    echo "[WARN] Unable to find matching R2: $R2; skipping $base"
    continue
  fi

  PREFIX="${OUTDIR}/${SAMPLE}"
  if [[ -f "${PREFIX}_S1.fastq.gz" && -f "${PREFIX}_S2.fastq.gz" ]]; then
    echo "[INFO] Output already exists, skipping: $SAMPLE"
    continue
  fi

  echo "[INFO] Deduplication: $SAMPLE"
  hts_SuperDeduper -1 "$R1" -2 "$R2" -f "$PREFIX" > "${PREFIX}.log" 2>&1

  if [[ -f "${PREFIX}_S1.fastq.gz" && -f "${PREFIX}_S2.fastq.gz" ]]; then
    ((pairs_done++))
  else
    echo "[ERROR] output missing：$SAMPLE"
  fi
done

echo "[INFO] completed sample pairs $pairs_done"

# --------- QC summary---------
SUMMARY="${OUTDIR}/dedup_summary.tsv"
echo -e "sample\traw_pairs(R1)\tdedup_pairs(R1)\tremoved_ratio" > "$SUMMARY"

for outR1 in "${OUTDIR}"/*_S1.fastq.gz; do
  [[ -e "$outR1" ]] || continue
  sample="$(basename "$outR1" | sed 's/_S1.fastq.gz//')"

  #  check R1 candiates 
  raw=""
  for cand in \
    "${INDIR}/${sample}_1_val_1.fq.gz"    "${INDIR}/${sample}_1_val_1.fastq.gz" \
    "${INDIR}/${sample}_1_val_1.fq"       "${INDIR}/${sample}_1_val_1.fastq" \
    "${INDIR}/${sample}_R1.fastq.gz"      "${INDIR}/${sample}_R1.fq.gz" \
    "${INDIR}/${sample}_1.fastq.gz"       "${INDIR}/${sample}_1.fq.gz" \
    "${INDIR}/${sample}_R1.fastq"         "${INDIR}/${sample}_R1.fq" \
    "${INDIR}/${sample}_1.fastq"          "${INDIR}/${sample}_1.fq" ; do
    [[ -f "$cand" ]] && raw="$cand" && break
  done
  [[ -z "$raw" ]] && { echo "[WARN] cannot find raw R1：$sample"; continue; }

  # calculate reads numbers
  if [[ "$raw" == *.gz ]]; then raw_n=$(zcat "$raw"   | awk 'END{print NR/4}'); else raw_n=$(awk 'END{print NR/4}' "$raw"); fi
  ded_n=$(zcat "$outR1" | awk 'END{print NR/4}')
  ratio=$(awk -v a="$raw_n" -v b="$ded_n" 'BEGIN{if(a>0) printf "%.4f",(a-b)/a; else print "NA"}')
  echo -e "${sample}\t${raw_n}\t${ded_n}\t${ratio}" >> "$SUMMARY"
done

echo "[INFO] summary completed:  $SUMMARY"
