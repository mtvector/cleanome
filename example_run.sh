#!/bin/bash
#SBATCH --job-name=cleanome_example
#SBATCH --output=/scratch/cleanome_genomes2/logs/%x.out
#SBATCH --error=/scratch/cleanome_genomes2/logs/%x.err
#SBATCH --time=48:00:00
#SBATCH --partition=celltypes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=48gb

set -euo pipefail

REPO_DIR="/code/cleanome"
JOB_DIR="/scratch/cleanome_genomes2"
GENOME_DIR="$JOB_DIR/ncbi_genomes"
ARC_OUT_DIR="$JOB_DIR/cellranger_arc"
RNA_OUT_DIR="$JOB_DIR/cellranger_rna"
ARC_SCRIPT_DIR="$JOB_DIR/submission_scripts"
RNA_SCRIPT_DIR="$JOB_DIR/rna_submission_scripts"
LOG_DIR="$JOB_DIR/logs"
STATS_CSV="$JOB_DIR/genome_info.csv"

SOURCE_GENOME_DIR="/data/genome_resources/cleanome_genomes2/ncbi_genomes"
SPECIES_LIST="$REPO_DIR/quick_list.txt"
PYTHON_BIN="${CLEANOME_PYTHON_BIN:-/opt/conda/bin/python}"
CLEANOME_SRC_DIR="${CLEANOME_SRC_DIR:-$REPO_DIR}"
CELLRANGER_ARC_BIN="/code/cellranger-arc-2.1.0/bin"
CELLRANGER_RNA_BIN="/code/cellranger-10.0.0/bin"

run_cleanome_python() {
    PYTHONPATH="$CLEANOME_SRC_DIR${PYTHONPATH:+:$PYTHONPATH}" "$PYTHON_BIN" "$@"
}

# stage_local_genomes() {
#     local src_dir dest_dir src base

#     mkdir -p "$GENOME_DIR"

#     for src_dir in "$SOURCE_GENOME_DIR"/*; do
#         [[ -d "$src_dir" ]] || continue

#         dest_dir="$GENOME_DIR/$(basename "$src_dir")"
#         mkdir -p "$dest_dir"

#         for src in "$src_dir"/*; do
#             [[ -f "$src" ]] || continue
#             base="$(basename "$src")"

#             case "$base" in
#                 *_splitchr*|*_debugged*|*_mito*)
#                     continue
#                     ;;
#             esac

#             ln -sfn "$src" "$dest_dir/$base"
#         done
#     done
# }

if [[ ! -x "$PYTHON_BIN" ]]; then
    echo "Unable to execute python at $PYTHON_BIN" >&2
    exit 1
fi

if [[ ! -d "$CLEANOME_SRC_DIR" ]]; then
    echo "Missing cleanome source tree at $CLEANOME_SRC_DIR" >&2
    exit 1
fi

if [[ ! -d "$SOURCE_GENOME_DIR" ]]; then
    echo "Missing source genome directory $SOURCE_GENOME_DIR" >&2
    exit 1
fi

if [[ ! -x "$CELLRANGER_ARC_BIN/cellranger-arc" ]]; then
    echo "Missing cellranger-arc at $CELLRANGER_ARC_BIN/cellranger-arc" >&2
    exit 1
fi

if [[ ! -x "$CELLRANGER_RNA_BIN/cellranger" ]]; then
    echo "Missing cellranger at $CELLRANGER_RNA_BIN/cellranger" >&2
    exit 1
fi

mkdir -p \
    "$JOB_DIR" \
    "$GENOME_DIR" \
    "$ARC_OUT_DIR" \
    "$RNA_OUT_DIR" \
    "$ARC_SCRIPT_DIR" \
    "$RNA_SCRIPT_DIR" \
    "$LOG_DIR"

# This creates a scratch-local mirror of the input genomes without copying the
# FASTA/GTF contents. Generated split/debugged files will then land on /scratch.
# stage_local_genomes

# If you want to download a fresh set of genomes instead of staging the local
# mirror above, comment out stage_local_genomes and uncomment this:
run_cleanome_python -m cleanome.download_genomes \
    --species_list "$SPECIES_LIST" \
    --genome_dir "$GENOME_DIR"

run_cleanome_python -m cleanome.get_genomes_and_stats \
    --genome_dir "$GENOME_DIR" \
    --stats_csv "$STATS_CSV"

run_cleanome_python -m cleanome.make_cellranger_arc_sh \
    --sh_scripts_dir "$ARC_SCRIPT_DIR" \
    --output_dir "$ARC_OUT_DIR" \
    --log_dir "$LOG_DIR" \
    --stats_csv "$STATS_CSV" \
    --cellranger_bin "$CELLRANGER_ARC_BIN" \
    --input_base_dir "$JOB_DIR" \
    --python_bin "$PYTHON_BIN" \
    --cleanome_src_dir "$CLEANOME_SRC_DIR"

run_cleanome_python -m cleanome.make_cellranger_rna_sh \
    --sh_scripts_dir "$RNA_SCRIPT_DIR" \
    --output_dir "$RNA_OUT_DIR" \
    --log_dir "$LOG_DIR" \
    --stats_csv "$STATS_CSV" \
    --cellranger_bin "$CELLRANGER_RNA_BIN"

cat <<EOF
Scratch-local cleanome example setup is complete.

Inputs staged under:
  $GENOME_DIR

Generated metadata:
  $STATS_CSV

Generated ARC scripts:
  $ARC_SCRIPT_DIR

Generated RNA scripts:
  $RNA_SCRIPT_DIR

Example species list for fresh downloads:
  $SPECIES_LIST
EOF
