"""
Author: Matthew Schmitz, Allen Institute, 2024

Downloads NCBI reference genomes and annotations for a list of species
usage:
python download_genomes.py {species_list_file} {genome_output_directory} {genome_info_csv}

Could also replace species list with a dataframe, or species dict as a dictionary of {taxid:scientific_name}
"""
import sys
import cleanome.genome_obtainment
from cleanome.genome_obtainment import *
import argparse

def main():
    parser = argparse.ArgumentParser(description="download genomes and generate summary stats.")
    
    parser.add_argument('-s', '--species_list', type=os.path.abspath, required=True, help='path to file of species names')
    parser.add_argument('-g', '--genome_dir', type=os.path.abspath, required=True, help='directory to download genomes to')
    parser.add_argument('-t','--taxid', action=argparse.BooleanOptionalAction, help='species list is taxids')
    args=parser.parse_args()
    
    genome_directory=args.genome_dir#'/home/matthew.schmitz/Matthew/genome/ncbi2'
    
    print('Downloading genomes...')
    with open(args.species_list, 'r') as file:
        species_list = [line.strip() for line in file]
    
    if args.taxid:
        species_dict = {species:species  for species in species_list}
    else:
        species_dict = {get_taxid(species):species  for species in species_list}
    
    for taxid in species_dict.keys():
        try:
            download_genome_data(str(taxid),genome_directory)
        except Exception as error:
            print(taxid," HAS FAILED", error)

if __name__ == "__main__":
    main()