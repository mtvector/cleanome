"""
Author: Matthew Schmitz, Allen Institute, 2024

Writes cellranger arc make scripts for all genomes in a table
usage:
Fill in variables and just run python make_cellranger_arc_sh.py
"""

import os
import re
import sys
import argparse
from cleanome._compat import block_broken_pyarrow

block_broken_pyarrow()

import pandas as pd


DEFAULT_INPUT_BASE_DIR = "/data/genome_resources/cleanome_genomes2"
DEFAULT_PYTHON_BIN = "/opt/conda/bin/python"
DEFAULT_CLEANOME_SRC_DIR = "/root/capsule/code/cleanome"
DEFAULT_CELLRANGER_ARC_BIN = "/code/cellranger-arc-2.1.0/bin"


def rebase_cleanome_path(path, input_base_dir):
    if pd.isna(path):
        return ""
    path = str(path).strip()
    if not path or path == "nan":
        return ""
    marker = f"{os.sep}cleanome_genomes2{os.sep}"
    if marker in path:
        return os.path.join(input_base_dir, path.split(marker, 1)[1])
    if os.path.isabs(path):
        return path
    return os.path.join(input_base_dir, path)


def strip_known_suffix(filename, pattern, fallback):
    filename = os.path.basename(filename) if filename else fallback
    return re.sub(pattern, "", filename)


def path_to_script_expr(path, base_dir, base_var_name):
    if not path:
        return ""
    abs_base_dir = os.path.abspath(base_dir)
    abs_path = os.path.abspath(path)
    if os.path.commonpath([abs_base_dir, abs_path]) == abs_base_dir:
        rel_path = os.path.relpath(abs_path, abs_base_dir)
        return f"${base_var_name}/{rel_path}"
    return abs_path

def main():

    parser = argparse.ArgumentParser(description="Make cellranger-arc reference generation bash scripts.")
    
    parser.add_argument('-s', '--sh_scripts_dir', type=os.path.abspath, required=True, help='Directory where the [Slurm] scripts will be saved (can just run them as normal bash)')
    parser.add_argument('-o', '--output_dir', type=os.path.abspath, required=True, help='Directory where cellranger-arc references will be output')
    parser.add_argument('-p', '--partition', type=str, default='celltypes',required=False, help='slurm partition')
    parser.add_argument('-c', '--cellranger_bin', type=os.path.abspath,default=DEFAULT_CELLRANGER_ARC_BIN, required=False, help='Path to cellranger-arc/bin')
    parser.add_argument('-l', '--log_dir', type=os.path.abspath, default=None, required=False, help='Path to log directory')
    parser.add_argument('-i', '--stats_csv', type=os.path.abspath, required=True, help='output csv file of genome metadata')
    parser.add_argument('-f', '--sif_dir', type=os.path.abspath, required=False, help='location of the singularity file to run mitofinder')
    parser.add_argument('--input_base_dir', type=os.path.abspath, default=DEFAULT_INPUT_BASE_DIR, required=False, help='Base directory for cleanome genome inputs')
    parser.add_argument('--python_bin', type=os.path.abspath, default=DEFAULT_PYTHON_BIN, required=False, help='Python interpreter used to run cleanome')
    parser.add_argument('--cleanome_src_dir', type=os.path.abspath, default=DEFAULT_CLEANOME_SRC_DIR, required=False, help='Cleanome source tree to place on PYTHONPATH')
    args=parser.parse_args()
    
    
    # Directory where the Slurm scripts will be saved
    sh_scripts_dir = args.sh_scripts_dir
    #arc output
    out_dir = args.output_dir
    #cellranger-arc/bin location
    cellranger_bin = args.cellranger_bin
    cellranger_root = os.path.dirname(cellranger_bin.rstrip(os.sep)) if cellranger_bin else ''
    log_dir = args.log_dir or os.path.join(out_dir, 'logs')
    # Ensure the scripts directory exists
    os.makedirs(sh_scripts_dir, exist_ok=True)
    
    # Load the CSV file which must have Species, Common Name, TaxID, FASTA Path and GTF Path columns
    df = pd.read_csv(args.stats_csv,sep=',',index_col=0)
    
    # Iterate over the rows of the DataFrame
    for index, row in df.iterrows():
        try:
            species = row['Species'].replace(' ', '_')
            common_name = row['Common Name'].replace(' ', '_')
            tax_spec=str(row['TaxID'])+"-"+row['Species']
            genome_assembly = tax_spec+"__"+re.sub(r'\.','-',re.sub(r'_genomic.+','',row['Genome Assembly Name']))
            print(genome_assembly)
            fasta_path = rebase_cleanome_path(row['FASTA Path'], args.input_base_dir)
            gtf_path = rebase_cleanome_path(row['GTF Path'], args.input_base_dir)
            source_fasta_expr = path_to_script_expr(fasta_path, args.input_base_dir, 'INPUT_BASE_DIR')
            source_gtf_expr = path_to_script_expr(gtf_path, args.input_base_dir, 'INPUT_BASE_DIR')
            split_fasta_filename = strip_known_suffix(fasta_path, r'(\.fasta|\.fna|\.fa)(\.gz)?$', 'genome') + '_splitchr.fa'
            split_gtf_filename = strip_known_suffix(gtf_path, r'(\.gtf|\.gff)(\.gz)?$', 'missing_input') + '_splitchr.gtf'
            gtf_debugged_filename = split_gtf_filename.replace('.gtf', '_debugged.gtf')
            config_file = f"{out_dir}/refseq_{species}.config"
            mito_id = str(row.get('Mitochondrial GenBank ID', '')).strip()
            mito_id = '' if mito_id == 'nan' else mito_id
            sif_dir = args.sif_dir if args.sif_dir else 'None'
        
            # Create the Slurm script content
            script_content = f"""#!/bin/bash
#SBATCH --job-name=mkref_{species}
#SBATCH --output={log_dir}/%x.out
#SBATCH --error={log_dir}/%x.err
#SBATCH --time=48:00:00
#SBATCH --partition={args.partition}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=48gb

INPUT_BASE_DIR="{args.input_base_dir}"
OUTPUT_BASE_DIR="{out_dir}"
NCBI_GENOMES_DIR="$INPUT_BASE_DIR/ncbi_genomes"
OUTPUT_NCBI_GENOMES_DIR="$OUTPUT_BASE_DIR/ncbi_genomes"
CELLRANGER_ARC_DIR="$OUTPUT_BASE_DIR"
LOG_DIR="$OUTPUT_BASE_DIR/logs"
CELLRANGER_ARC_BIN_DIR="{cellranger_bin}"
CELLRANGER_ARC_ROOT_DIR="{cellranger_root}"
GTF_TO_GENE_INDEX_BIN="${{CELLRANGER_ARC_ROOT_DIR}}/lib/bin/gtf_to_gene_index"
PYTHON_BIN="${{CLEANOME_PYTHON_BIN:-{args.python_bin}}}"
CLEANOME_SRC_DIR="${{CLEANOME_SRC_DIR:-{args.cleanome_src_dir}}}"

SOURCE_FASTA="{source_fasta_expr}"
SOURCE_GTF="{source_gtf_expr}"
SPECIES_DIR="$OUTPUT_NCBI_GENOMES_DIR/{tax_spec}"
SPLIT_FASTA="$SPECIES_DIR/{split_fasta_filename}"
SPLIT_GTF="$SPECIES_DIR/{split_gtf_filename}"
DEBUGGED_GTF="$SPECIES_DIR/{gtf_debugged_filename}"
CONFIG_PATH="$CELLRANGER_ARC_DIR/refseq_{species}.config"
REF_NAME="{genome_assembly}"
SIF_DIR="{sif_dir}"

mkdir -p "$LOG_DIR" "$CELLRANGER_ARC_DIR" "$SPECIES_DIR"

run_cleanome_python() {{
    PYTHONPATH="$CLEANOME_SRC_DIR${{PYTHONPATH:+:$PYTHONPATH}}" "$PYTHON_BIN" "$@"
}}

if [[ ! -x "$PYTHON_BIN" ]]; then
    echo "Unable to execute python at $PYTHON_BIN" >&2
    exit 1
fi
if [[ ! -d "$CLEANOME_SRC_DIR" ]]; then
    echo "Missing cleanome source tree at $CLEANOME_SRC_DIR" >&2
    exit 1
fi
if [[ ! -x "$CELLRANGER_ARC_BIN_DIR/cellranger-arc" ]]; then
    echo "Missing cellranger-arc at $CELLRANGER_ARC_BIN_DIR/cellranger-arc" >&2
    exit 1
fi
if [[ ! -x "$GTF_TO_GENE_INDEX_BIN" ]]; then
    echo "Missing gtf_to_gene_index at $GTF_TO_GENE_INDEX_BIN" >&2
    exit 1
fi
if [[ -z "$SOURCE_FASTA" || ! -f "$SOURCE_FASTA" ]]; then
    echo "Missing source FASTA $SOURCE_FASTA" >&2
    exit 1
fi
if [[ -z "$SOURCE_GTF" || ! -f "$SOURCE_GTF" ]]; then
    echo "Missing source GTF $SOURCE_GTF" >&2
    exit 1
fi

run_cleanome_python - <<EOF
from cleanome.chromosome_splitter import ChromosomeSplitter
ChromosomeSplitter(
    "$SOURCE_FASTA",
    "$SOURCE_GTF",
    "$SPLIT_FASTA",
    "$SPLIT_GTF",
    length_threshold=5e8
)
EOF

# — prepare split files —
TARGET_FASTA="$SPLIT_FASTA"
TARGET_GTF="$SPLIT_GTF"

# — if the user provided a mito accession, run add_mito now! —
mito_id="{mito_id}"
if [[ -n "$mito_id" && "$mito_id" != "nan" ]]; then
    ref_gbk="$(dirname "$SOURCE_FASTA")/$mito_id.fa"
    if [[ -f "$ref_gbk" ]]; then
        run_cleanome_python -m cleanome.add_mito \\
            -a "$mito_id" \\
            -r "$ref_gbk" \\
            -g "$TARGET_FASTA" \\
            -t "$TARGET_GTF" \\
            -s "$SIF_DIR" \\
            -n "{species}"

        # now point at the mito-augmented files
        TARGET_FASTA="${{TARGET_FASTA%.*}}_mito.fasta"
        TARGET_GTF="${{TARGET_GTF%.gtf}}_mito.gtf"
    else
        echo "Warning: missing GBK $ref_gbk; skipping add_mito"
    fi
fi

# — Step 1 (new): Debug _that_ GTF so the mito lines are clean too —
DEBUGGED_GTF="${{TARGET_GTF%.gtf}}_debugged.gtf"
run_cleanome_python -m cleanome.debug_gtf "$TARGET_GTF" "$DEBUGGED_GTF"
FINAL_FASTA="$TARGET_FASTA"
FINAL_GTF="$DEBUGGED_GTF"

cd "$CELLRANGER_ARC_DIR"

# Step 2: Create the config file for cellranger-arc
cat > "$CONFIG_PATH" <<EOF
{{
    organism: "{common_name}"
    genome: ["{genome_assembly}"]
    input_fasta: ["$FINAL_FASTA"]
    input_gtf: ["$FINAL_GTF"]
}}
EOF


# Step 3: Run cellranger-arc mkref
"$CELLRANGER_ARC_BIN_DIR/cellranger-arc" mkref --config="$CONFIG_PATH"
cp "$CONFIG_PATH" "$REF_NAME"

"$GTF_TO_GENE_INDEX_BIN" "$CELLRANGER_ARC_DIR/$REF_NAME" test.json
            """
            
            # Save the script to a file
            script_filename = os.path.join(sh_scripts_dir, f"slurm_{species}.sh")
            with open(script_filename, 'w') as script_file:
                script_file.write(script_content)
            os.chmod(script_filename, 0o755)
        
            print(f"Slurm script for {species} saved to {script_filename}")
        except Exception as exc:
            print(f"Failed to generate script for {row.get('Species', index)}: {exc}", file=sys.stderr)
    print("All Slurm scripts have been generated.")
    print('submit with: for script in *.sh; do sbatch "$script"; done ')

if __name__ == "__main__":
    main()
