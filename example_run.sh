#!/bin/bash
#SBATCH --job-name=mkrefs
#SBATCH --output=/home/matthew.schmitz/log/testrefs.out
#SBATCH --error=/home/matthew.schmitz/log/testrefs.err
#SBATCH --time=48:00:00
#SBATCH --partition=celltypes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=48gb

source ~/.bashrc
conda activate cleanome

CLEANOME_DIR=~/utils/cleanome
JOB_DIR=/home/matthew.schmitz/Matthew/genome/cleanome_test
mkdir -p $JOB_DIR

cd $JOB_DIR
mkdir -p ncbi_genomes
mkdir -p cellranger_arc
mkdir -p submission_scripts
mkdir -p logs

download_genomes --species_list ~/species_list.txt \
--genome_dir $JOB_DIR/ncbi_genomes/ 

get_genomes_and_stats --genome_dir $JOB_DIR/ncbi_genomes/ \
--stats_csv $JOB_DIR/genome_info.csv -c

make_cellranger_arc_sh --sh_scripts_dir $JOB_DIR/submission_scripts/ \
--output_dir $JOB_DIR/cellranger_arc \
--log_dir $JOB_DIR/logs \
--stats_csv $JOB_DIR/genome_info.csv \
--cellranger_bin /allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/utils/cellranger-arc-2.0.2/bin
