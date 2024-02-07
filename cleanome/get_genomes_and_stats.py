"""
Author: Matthew Schmitz, Allen Institute, 2024

Downloads NCBI reference genomes and annotations for a list of species
usage:
python get_genomes_and_stats.py {species_list_file} {genome_output_directory} {genome_info_csv}

Could also replace species list with a dataframe, or species dict as a dictionary of {taxid:scientific_name}
"""
import sys
import cleanome.genome_obtainment
from cleanome.genome_obtainment import *
import argparse

def main():
    parser = argparse.ArgumentParser(description="download genomes and generate summary stats.")
    
    parser.add_argument('-g', '--genome_dir', type=os.path.abspath, required=True, help='directory to download genomes to')
    parser.add_argument('-o', '--stats_csv', type=os.path.abspath, required=True, help='output csv file of genome metadata')
    parser.add_argument('-c','--calculate_stats', action=argparse.BooleanOptionalAction, help='Calculate transcriptome stats from GTF? (can take a long time)')

    args=parser.parse_args()
    
    genome_directory=args.genome_dir#'/home/matthew.schmitz/Matthew/genome/ncbi2'
    genome_csv=args.stats_csv#'./AllNCBIGenomesData.csv'
    
    print('Collecting Genome Info...')
    ncbi_data = list(process_directory(genome_directory,gtf=args.calculate_stats))
    
    ncbi_data = [x for x in ncbi_data if x is not None]
    
    # Create DataFrames
    ncbi_df = pd.DataFrame(ncbi_data)
    
    combined_df.to_csv(ncbi_df)

if __name__ == "__main__":
    main()
