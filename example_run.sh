#!/bin/bash
#SBATCH --job-name=mkrefs
#SBATCH --output=/home/matthew.schmitz/log/testrefs.out
#SBATCH --error=/home/matthew.schmitz/log/testrefs.out
#SBATCH --time=48:00:00
#SBATCH --partition=celltypes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=48gb

PYTHON_BIN="${CLEANOME_PYTHON_BIN:-$(command -v python)}"
if [[ -z "$PYTHON_BIN" ]]; then
    echo "Unable to locate python; set CLEANOME_PYTHON_BIN" >&2
    exit 1
fi

JOB_DIR=/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/genomes/cleanome_genomes2
mkdir -p $JOB_DIR

cd $JOB_DIR
mkdir -p ncbi_genomes
mkdir -p cellranger_arc
mkdir -p cellranger_rna
mkdir -p submission_scripts
mkdir -p rna_submission_scripts
mkdir -p logs

# download_genomes --species_list ~/utils/cleanome/species_list.txt \
# --genome_dir $JOB_DIR/ncbi_genomes/ 

"$PYTHON_BIN" -m cleanome.get_genomes_and_stats --genome_dir $JOB_DIR/ncbi_genomes/ \
--stats_csv $JOB_DIR/genome_info.csv # -c

# Here you would locate mitochondria genomes and add the to genome_info.csv if you wish
# Or add them independently with add_mito

"$PYTHON_BIN" -m cleanome.make_cellranger_arc_sh --sh_scripts_dir $JOB_DIR/submission_scripts/ \
--output_dir $JOB_DIR/cellranger_arc \
--log_dir $JOB_DIR/logs \
--stats_csv $JOB_DIR/genome_info.csv \
--cellranger_bin /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/cellranger-arc-2.1.0/bin

"$PYTHON_BIN" -m cleanome.make_cellranger_rna_sh --sh_scripts_dir $JOB_DIR/rna_submission_scripts/ \
--output_dir $JOB_DIR/cellranger_rna \
--log_dir $JOB_DIR/logs \
--stats_csv $JOB_DIR/genome_info.csv \
--cellranger_bin /allen/programs/celltypes/workgroups/rnaseqanalysis/bicore/tools/cellranger-8.0.0/bin
