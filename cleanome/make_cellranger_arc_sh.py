"""
Author: Matthew Schmitz, Allen Institute, 2024

Writes cellranger arc make scripts for all genomes in a table
usage:
Fill in variables and just run python make_cellranger_arc_sh.py
"""

import pandas as pd
import os
import re
import sys
from cleanome.chromosome_splitter import *
import argparse

def main():

    parser = argparse.ArgumentParser(description="Make cellranger-arc reference generation bash scripts.")
    
    parser.add_argument('-s', '--sh_scripts_dir', type=os.path.abspath, required=True, help='Directory where the [Slurm] scripts will be saved (can just run them as normal bash)')
    parser.add_argument('-o', '--output_dir', type=os.path.abspath, required=True, help='Directory where cellranger-arc references will be output')
    parser.add_argument('-p', '--partition', type=str, default='celltypes',required=False, help='slurm partition')
    parser.add_argument('-c', '--cellranger_bin', type=os.path.abspath,default='', required=False, help='Path to cellranger-arc/bin ($PATH finding is unreliable, but you can rely on it if you want)')
    parser.add_argument('-l', '--log_dir', type=os.path.abspath ,default='./', required=False, help='Path to cellranger-arc/bin ($PATH finding is unreliable, but you can rely on it if you want)')
    parser.add_argument('-i', '--stats_csv', type=os.path.abspath, required=True, help='output csv file of genome metadata')
    parser.add_argument('-f', '--sif_dir', type=os.path.abspath, required=False, help='location of the singularity file to run mitofinder')
    args=parser.parse_args()
    
    
    # Directory where the Slurm scripts will be saved
    sh_scripts_dir = args.sh_scripts_dir
    #arc output
    out_dir = args.output_dir
    #cellranger-arc/bin location
    cellranger_bin = args.cellranger_bin
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
            genome_assembly = tax_spec+"__"+re.sub('\.','-',re.sub('_genomic.+','',row['Genome Assembly Name']))
            print(genome_assembly)
            fasta_path = row['FASTA Path']
            gtf_path = row['GTF Path']
            split_fasta_path=re.sub('\.fasta|\.fna|\.fa|\.gz','',fasta_path)+'_splitchr.fa'
            split_gtf_path=re.sub('\.gtf|\.gff|\.gz','',gtf_path)+'_splitchr.gtf'
            gtf_debugged_filename = os.path.basename(split_gtf_path).replace('.gtf', '_debugged.gtf')
            gtf_debugged_path = os.path.join(os.path.dirname(split_gtf_path), gtf_debugged_filename)
            config_file = f"{out_dir}/refseq_{species}.config"
        
            # Create the Slurm script content
            script_content = f"""#!/bin/bash
#SBATCH --job-name=mkref_{species}
#SBATCH --output={args.log_dir}/mkref_{species}.out
#SBATCH --error={args.log_dir}/mkref_{species}.err
#SBATCH --time=48:00:00
#SBATCH --partition={args.partition}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=48gb

# Load Python environment 
source ~/.bashrc
source ~/.zshrc

conda activate cleanome

python - <<EOF
from cleanome.chromosome_splitter import ChromosomeSplitter
ChromosomeSplitter(
    "{fasta_path}",
    "{gtf_path}",
    "{split_fasta_path}",
    "{split_gtf_path}",
    length_threshold=5e8
)
EOF

# — prepare split files —
TARGET_FASTA="{split_fasta_path}"
TARGET_GTF="{split_gtf_path}"
SIF_DIR="{args.sif_dir}"

# — if the user provided a mito accession, run add_mito now! —
mito_id="{row.get('Mitochondrial GenBank ID','')}"
if [[ -n "$mito_id" && "$mito_id" != "nan" ]]; then
    ref_gbk="$(dirname "{fasta_path}")/$mito_id.fa"
    if [[ -f "$ref_gbk" ]]; then
        python -m cleanome.add_mito \\
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
debug_gtf "$TARGET_GTF" "$DEBUGGED_GTF"
FINAL_FASTA="$TARGET_FASTA"
FINAL_GTF="$DEBUGGED_GTF"

cd {out_dir}

# Step 2: Create the config file for cellranger-arc
echo "{{
    organism: \\"{common_name}\\"
    genome: [\\"{genome_assembly}\\"]
    input_fasta: [\\"$FINAL_FASTA\\"]
    input_gtf: [\\"$FINAL_GTF\\"]
}}" > {config_file}


# Step 3: Run cellranger-arc mkref
export PATH=$PATH:{cellranger_bin}
cellranger-arc mkref --config={config_file}
cp {config_file} {genome_assembly}

{cellranger_bin}/gtf_to_gene_index {out_dir}/{genome_assembly} test.json 
            """
            
            # Save the script to a file
            script_filename = os.path.join(sh_scripts_dir, f"slurm_{species}.sh")
            with open(script_filename, 'w') as script_file:
                script_file.write(script_content)
        
            print(f"Slurm script for {species} saved to {script_filename}")
        except:
            pass
    print("All Slurm scripts have been generated.")
    print('submit with: for script in *.sh; do sbatch "$script"; done ')

if __name__ == "__main__":
    main()
