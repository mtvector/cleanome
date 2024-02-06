"""
Author: Matthew Schmitz, Allen Institute, 2024

Downloads NCBI reference genomes and annotations for a list of species
usage:
python download_genomes.py {species_list_file} {genome_output_directory} {genome_info_csv}

Could also replace species list with a dataframe, or species dict as a dictionary of {taxid:scientific_name}
"""
import sys
sys.path.append('/home/matthew.schmitz/utils/cleanome/cleanome')
import genome_obtainment
from genome_obtainment import *

import argparse

parser = argparse.ArgumentParser(description="download genomes and generate summary stats.")

parser.add_argument('-s', '--species_list', type=os.path.abspath, required=True, help='path to file of species names')
parser.add_argument('-g', '--genome_dir', type=os.path.abspath, required=True, help='directory to download genomes to')
parser.add_argument('-o', '--stats_csv', type=os.path.abspath, required=True, help='output csv file of genome metadata')
parser.add_argument('-d', '--download', type=bool, default=True, required=False, help='Download genomes (false just makes csv for output dir)')
parser.add_argument('-t', '--taxid', type=bool, default=False, required=False, help='Species list is taxids?')
parser.add_argument('-c', '--calculate_stats', type=bool, default=False, required=False, help='Calculate transcriptome stats from GTF? (can take a long time)')

genome_directory=parser.genome_dir#'/home/matthew.schmitz/Matthew/genome/ncbi2'
genome_csv=parser.stats_csv#'./AllNCBIGenomesData.csv'

if parser.download:
    with open(parser.species_list, 'r') as file:
        species_list = [line.strip() for line in file]
    
    if parser.tax_id:
        species_dict = {species:species  for species in species_list}
    else:
        species_dict = {get_taxid(species):species  for species in species_list}
    
    for taxid in species_dict.keys():
        download_genome_data(str(taxid),genome_directory)

ncbi_data = list(process_directory(genome_directory,gtf=parser.calculate_stats))

ncbi_data = [x for x in ncbi_data if x is not None]

# Create DataFrames
ncbi_df = pd.DataFrame(ncbi_data)

# Combine DataFrames
combined_df = pd.concat([ ncbi_df, other_df])
print(combined_df)
combined_df.to_csv(genome_csv)