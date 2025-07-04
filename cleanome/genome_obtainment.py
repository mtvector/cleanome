import sys
import subprocess
import re
import shutil
import requests
import os
import re
from ftplib import FTP
import pandas as pd
from Bio import SeqIO
from Bio import Entrez
import gtfparse
import gzip
from ete3 import Tree
import copy
import hashlib


def compute_md5(path, block_size=65536):
    """
    Compute the MD5 checksum of a file on disk.
    """
    md5 = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(block_size), b''):
            md5.update(chunk)
    return md5.hexdigest()

def max_column_match(v,df):
    max_count=0
    max_ind=''
    for c in df.columns:
        matched_count=df[c].isin(v).sum()
        if matched_count>max_count:
            max_count=matched_count
            max_ind=c
    return(max_ind)

def generate_chromsizes(fasta_file, output_file):
    with open(fasta_file, 'r') as infile, open(output_file, 'w') as outfile:
        for record in SeqIO.parse(infile, 'fasta'):
            sequence_name = record.id
            sequence_length = len(record.seq)
            outfile.write(f'{sequence_name}\t{sequence_length}\n')


def get_taxid(species_name):
    """Get the taxid for a species name using NCBI's Entrez Utilities."""
    Entrez.email = "your-email@example.com"  # Provide your email address
    handle = Entrez.esearch(db="taxonomy", term=species_name)
    record = Entrez.read(handle)
    handle.close()
    
    # Return the first TaxID found
    if record["IdList"]:
        return record["IdList"][0]
    else:
        return None


def get_scientific_name(taxid):
    """Get the scientific name for a given TaxID using NCBI's Entrez Utilities."""
    Entrez.email = "your-email@example.com"  # Provide your email address
    handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    
    # Return the scientific name found
    if records:
        return records[0]["ScientificName"]
    else:
        return None


#######################NCBI REFERENCE PROCESSING########################


def sanitize_filename(filename):
    # Replace spaces with underscores
    filename = filename.replace(" ", "_")

    # Replace special characters with dashes
    # Define a regex pattern for unwanted characters
    pattern = r'[^\w\-_.]'
    filename = re.sub(pattern, '-', filename)

    return filename

def run_ncbi_datasets_command(command):
    """Run an NCBI datasets CLI command and handle output."""
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running command {' '.join(command)}: {e}")

def download_genome_data(species_taxid, target_directory):
    """Download genome data using NCBI datasets CLI."""
    # Create the species-specific directory
    handle = Entrez.esummary(db="taxonomy", id=species_taxid)
    record = Entrez.read(handle)
    scientific_name = record[0]['ScientificName']
    taxid_species=str(species_taxid)+"-"+sanitize_filename(scientific_name)
    species_dir = os.path.join(target_directory, taxid_species)
    if os.path.exists(species_dir):
        print('exists')
        return None
    os.makedirs(species_dir, exist_ok=True)
    zip_file=os.path.join(target_directory,f"{species_taxid}.zip")
    # Download the genome data
    command = [
        "datasets", 
        "download", 
        "genome", 
        "taxon", 
        species_taxid, 
        "--reference", 
        "--assembly-source",
        "RefSeq",
        "--include", 
        "genome,gtf,seq-report",
        "--filename", 
        zip_file
    ]
    run_ncbi_datasets_command(command)

    # Unzip the downloaded file
    subprocess.run(["unzip", zip_file, "-d", species_dir], stdout=open(os.devnull, 'wb'))

    # Move files from their subdirectories to the main species directory
    for root, dirs, files in os.walk(species_dir):
        for file in files:
            file_path = os.path.join(root, file)
            if file_path.endswith(('.fna','.fa','.fasta', '.gtf','.json','.gff','.gff3','.jsonl')):
                shutil.move(file_path, species_dir)



##############################ENSEMBL 2023 Processing################################################



def sanitize_filename(filename):
    # Replace spaces with underscores and special characters with dashes
    filename = filename.replace(" ", "_")
    pattern = r'[^\w\-_.]'
    return re.sub(pattern, '-', filename)

def download_ftp_file(ftp_url, destination_filename):
    # Parse and download from FTP URL
    ftp_host, ftp_path = ftp_url.replace("ftp://", "").split("/", 1)
    with FTP(ftp_host) as ftp:
        ftp.login()
        os.makedirs(os.path.dirname(destination_filename), exist_ok=True)
        try:
            with open(destination_filename, 'wb') as local_file:
                ftp.retrbinary(f'RETR {ftp_path}', local_file.write)
        except Exception as e:
            print(f"Error downloading {ftp_url}: {e}")

def get_species_name(tax_id):
    """ Fetch species name using taxonomic ID from Ensembl """
    url = f"https://rest.ensembl.org/taxonomy/id/{tax_id}?content-type=application/json"
    response = requests.get(url)
    return response.json()['scientific_name']

def find_matching_files(ftp, path, patterns):
    """ Find files matching given patterns """
    ftp.cwd(path)
    files = ftp.nlst()
    matched_files = [file for file in files if any(pattern in file for pattern in patterns)]
    return matched_files

def download_ensembl_files(tax_id, path,redownload=False):
    species_name = get_species_name(tax_id)
    species_path = species_name.lower().replace(' ', '_')

    # Fetch genome information from Ensembl
    url = f"https://rest.ensembl.org/info/genomes/{species_path}?content-type=application/json"
    print(url)
    response = requests.get(url)
    data = response.json()
    assembly_name = data['assembly_name']

    # Prepare directory name
    dir_name = sanitize_filename(f"{tax_id}-{species_name}")
    dir_path = os.path.join(path, dir_name)
    os.makedirs(dir_path, exist_ok=True)

    # Prepare FTP paths
    ensembl_ftp_base = "ftp.ensembl.org"
    fasta_ftp_dir = f"/pub/release-110/fasta/{species_path}/dna/"
    gtf_ftp_dir = f"/pub/release-110/gtf/{species_path}/"

    with FTP(ensembl_ftp_base) as ftp:
        ftp.login()

        # Find matching FASTA and GTF files
        fa_patterns = [".dna_sm.toplevel.fa.gz"]
        gtf_patterns = [".gtf.gz"]
        unwanted_patterns = ["chr", "hapl", "abinitio"]

        fasta_files = find_matching_files(ftp, fasta_ftp_dir, fa_patterns)
        gtf_files = find_matching_files(ftp, gtf_ftp_dir, gtf_patterns)
        fasta_files = [f for f in fasta_files if not any(p in f for p in unwanted_patterns)]
        gtf_files = [f for f in gtf_files if not any(p in f for p in unwanted_patterns)]
        if not redownload and os.path.exists(os.path.join(dir_path, fasta_files[0])):
            return None
        # Download files
        if fasta_files:
            download_ftp_file(f"ftp://{ensembl_ftp_base}{fasta_ftp_dir}{fasta_files[0]}", os.path.join(dir_path, fasta_files[0]))
        if gtf_files:
            download_ftp_file(f"ftp://{ensembl_ftp_base}{gtf_ftp_dir}{gtf_files[0]}", os.path.join(dir_path, gtf_files[0]))


############################ASSESSING REFERENCE QUALITY##################################

import os
import json

def find_file_accession_and_assembly_name(directory, filename='assembly_data_report.jsonl'):
    # Walk through all directories and files in the given directory
    for root, dirs, files in os.walk(directory):
        if filename in files:
            # Construct the full path to the file
            file_path = os.path.join(root, filename)
            try:
                # Open and read the JSONL file
                with open(file_path, 'r') as file:
                    for line in file:
                        # Parse each line (assuming each line is a valid JSON object)
                        data = json.loads(line)
                        # Extract the accession value
                        accession = data.get('accession', 'No accession found')
                        # Extract the assemblyName value
                        assembly_name = data.get('assemblyInfo', {}).get('assemblyName', 'No assembly name found')
                        return accession, assembly_name, file_path
            except Exception as e:
                print(f"Error reading or parsing {file_path}: {e}")
                return None, None, None
    # Return None if the file is not found
    return None, None, None


def is_gzipped(file_path):
    """ Check if the file is gzipped """
    with open(file_path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def calculate_N50(fasta_path):
    file_open_func = gzip.open if is_gzipped(fasta_path) else open
    with file_open_func(fasta_path, "rt") as handle:
        lengths = [len(rec) for rec in SeqIO.parse(handle, "fasta")]
    total_len = sum(lengths)
    lengths.sort(reverse=True)
    cumsum = 0
    for length in lengths:
        cumsum += length
        if cumsum >= total_len / 2:
            return length
    return 0

"""
def count_genes_transcripts(gtf_path):
    db = gffutils.create_db(gtf_path, dbfn=":memory:", force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    num_genes = len(set(db.features_of_type('gene_id')))
    num_transcripts = len(set(db.features_of_type('transcript_id')))
    return num_genes, num_transcripts
"""

def count_genes_transcripts(gtf_path):
    try:
        gtf_df = gtfparse.read_gtf(gtf_path)
        num_genes = len(gtf_df['gene_id'].unique())
        try:
            num_transcripts = len(gtf_df['transcript_id'].unique())
        except:
            num_transcripts = num_genes
        return num_genes, num_transcripts
    except:
        return 0,0


EXCLUDE_SUFFIXES = ('_debugged', '_splitchr','_mito')

def find_original_file(directory, pattern):
    """
    Return the first filename in `directory` matching `pattern` 
    but not containing any of EXCLUDE_SUFFIXES.
    """
    for fn in os.listdir(directory):
        if re.search(pattern, fn) and not any(s in fn for s in EXCLUDE_SUFFIXES):
            return fn
    return None

def find_files_in_directory(directory, file_extension):
    for file in os.listdir(directory):
        if re.search(file_extension, file):
            return file
    return None


def build_ncbi_ftp_link(accession, assembly_name, ext):
    """
    accession: e.g. 'GCF_030014405.1'
    assembly_name: e.g. 'ASM3001440v1'
    ext: one of 'genomic.fna.gz' or 'genomic.gtf.gz'
    """
    acc_no, version = accession.split('.')
    prefix, num = acc_no.split('_')
    parts = [num[i:i+3] for i in range(0, len(num), 3)]
    dir_path = "/".join([prefix] + parts + [f"{acc_no}.{version}_{assembly_name}"])
    fname = f"{acc_no}.{version}_{assembly_name}_{ext}"
    return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{dir_path}/{fname}"


def process_directory(dir_path, gtf_stats=False):
    for root, dirs, files in os.walk(dir_path):
        for d in dirs:
            dir_full = os.path.join(root, d)
            m = re.search(r'([0-9]+-[A-Za-z_]+)', dir_full)
            if not m: 
                continue
            taxid, species_tag = m.group(0).split('-',1)
            # find only original fasta/gtf
            fasta_file = find_original_file(dir_full, r'\.f(n)?a(\.gz)?$')
            if not fasta_file:
                continue
            fasta_path = os.path.join(dir_full, fasta_file)
            accession, assembly_name, _ = find_file_accession_and_assembly_name(dir_full)
            assembly_name = assembly_name or fasta_file
            n50 = calculate_N50(fasta_path)
            gtf_file = find_original_file(dir_full, r'\.gtf(\.gz)?$')
            if gtf_file:
                gtf_path = os.path.join(dir_full, gtf_file)
            else:
                gtf_path = None
            fasta_md5 = compute_md5(fasta_path) if fasta_file else None
            gtf_md5   = compute_md5(gtf_path)   if gtf_file   else None
            if gtf_stats and gtf_file:
                num_genes, num_transcripts = count_genes_transcripts(gtf_path)            
            else:
                num_genes = num_transcripts = None
            
            stat = {
                'Species': species_tag,
                'TaxID': taxid,
                'Common Name': species_tag.replace('_',' '),
                'Genome Assembly Name': assembly_name,
                'Genome Accession': accession,
                'N50': n50,
                'Number of Genes': num_genes,
                'Number of Transcripts': num_transcripts,
                'FASTA File': fasta_file,
                'FASTA Path': fasta_path,
                'FASTA MD5': fasta_md5,
                'GTF File': gtf_file,
                'GTF Path': gtf_path,
                'GTF MD5': gtf_md5,
                'Mitochondrial GenBank ID': '',
                'FTP FASTA URL': build_ncbi_ftp_link(accession, assembly_name, 'genomic.fna.gz') if accession else '',
                'FTP GTF URL':   build_ncbi_ftp_link(accession, assembly_name, 'genomic.gtf.gz')   if accession else '',
            }
            yield stat

def get_genera_from_species_list(species_list):
    """Extract genera from a list of species' scientific names."""
    return set(species.split()[0] for species in species_list)

def prune_tree(newick_tree_path, genera):
    """Prune the tree to retain only species from the specified genera."""
    tree = Tree(newick_tree_path,format=1)
    original_tree = copy.deepcopy(tree)
    for leaf in tree.iter_leaves():
        genus = leaf.name.split()[0]
        if genus not in genera:
            leaf.delete()
    return tree,original_tree

def get_species_not_in_tree(species_list, tree):
    # Extract species names from the tree
    species_in_tree = set(leaf.name for leaf in tree.iter_leaves())

    # Find species not in the tree
    species_not_in_tree = [species for species in species_list if species not in species_in_tree]

    return species_not_in_tree

def add_node_to_common_ancestor(tree, new_node_name, existing_node_names):
    all_nodes = [node.name for node in tree.traverse()]
    if new_node_name in all_nodes:
        print(new_node_name,'is in the tree already!')
        return None
    # Find the common ancestor of the existing nodes
    common_ancestor = tree.get_common_ancestor(existing_node_names)

    # Calculate the mean branch length to the common ancestor
    branch_lengths = [tree.get_distance(name,common_ancestor) for name in existing_node_names]
    mean_branch_length = sum(branch_lengths) / len(branch_lengths)

    # Create the new node and add it to the common ancestor
    new_node = common_ancestor.add_child(name=new_node_name, dist=mean_branch_length)
    
    return new_node

def filter_genome_dataframe(df):
    # Sort the dataframe first by 'Species', then by 'Number of Transcripts' (descending), and then by 'N50' (descending)
    df_sorted = df.sort_values(by=['Species', 'Number of Transcripts', 'N50'], ascending=[True, False, False])

    # Drop duplicates, keeping the first occurrence (which is the one with the highest 'Number of Transcripts' or 'N50' in case of a tie)
    df_filtered = df_sorted.drop_duplicates(subset='Species', keep='first')

    return df_filtered

def write_tree_and_dataframe_to_file(tree, df, filename, col1, col2):
    """
    Writes a Newick tree, a space, and two DataFrame columns to a file.

    :param tree: Newick tree string.
    :param df: Pandas DataFrame.
    :param filename: Output file name.
    :param col1: First column name from DataFrame.
    :param col2: Second column name from DataFrame.
    """
    with open(filename, 'w') as f:
        # Write the Newick tree
        f.write(tree + '\n\n')  # Newick tree followed by a newline and a space

        # Write the specified DataFrame columns
        df[[col1, col2]].to_csv(f, sep='\t', index=False,header=None)


###################################Ortho Tables#################


def download_biomart_table(dataset, filename):
    #all_attributes = [x for x in dataset.attributes.keys() if 'hsapiens' in x]
    all_attributes = [
        'ensembl_gene_id',
        'external_gene_name',
        'hsapiens_homolog_associated_gene_name',
        'hsapiens_homolog_ensembl_gene',
        'hsapiens_homolog_orthology_confidence'
    ]

    response = dataset.search({
        'attributes': all_attributes
    }, header=1)  # header=1 will include the column names

    with open(filename, 'wb') as f:
        f.write('\t'.join(all_attributes).encode('ascii')+b'\n')
        for line in response.iter_lines():
            f.write(line + b'\n')

def find_matching_rows(strings, df,colnames=None):
    matches = {}
    if colnames is None:
        colnames=df.columns
    use_df=df.loc[:,colnames]
    for s in tqdm.tqdm(strings):
        mask = use_df.isin([s]).any(axis=1)
        matches[s] = df[mask]
    return matches

def get_human_orthologs(species_name, identifiers, cache_path):
    #Returns a dictionary of {original_key: dataframe of matched rows}
    server = BiomartServer("http://www.ensembl.org/biomart")
    dataset_name = f"{species_name.lower()}_gene_ensembl"
    
    # Construct cache_path based on dataset_name
    cache_name = os.path.join(cache_path,dataset_name + "_table.txt")
    
    dataset = server.datasets[dataset_name]

    # Check cache and download if necessary
    if not os.path.exists(cache_name):
        print("Cache not found. Downloading BioMart table...")
        download_biomart_table(dataset, cache_name)
        print(f"Downloaded and saved to {cache_name}")

    # Read the table into a DataFrame
    df = pd.read_csv(cache_name, sep='\t', dtype=str)
    
    # Filter using Ensembl IDs and gene symbols and add 'original_identifier' column
    #filtered_df = df[df['ensembl_gene_id'].isin(identifiers) | df['external_gene_name'].isin(identifiers)].copy()
    #filtered_df['original_identifier'] = filtered_df.apply(lambda row: row['ensembl_gene_id'] if row['ensembl_gene_id'] in identifiers else row['external_gene_name'], axis=1)
    #orthologs = {row['original_identifier']: row.drop('original_identifier').to_dict() for _, row in filtered_df.iterrows()}

    orthologs=find_matching_rows(identifiers, df,df.columns[df.columns.str.contains('gene')])        
    return orthologs
