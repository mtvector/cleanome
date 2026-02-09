# cleanome

**Cleanome** is a Python-based helper toolkit designed to standardize and clean genome annotation files (especially GTFs) so they can be consumed by genomics reference builders like Cellranger (and Cellranger-ARC). At a high level, running a genome through Cleanome involves the following modifications and preprocessing steps: 
*	**Downloading core files**: It fetches the FASTA sequence, GTF annotation, and associated assembly metadata for a list of species you provide (by scientific name or taxonomy ID), centralizing them in a working directory. 
*	**Diagnostic statistics**: It computes basic assembly and annotation statistics across all genomes (e.g., contig sizes, number of genes/transcripts), producing a summary table used later for conditional corrections. 
* **GTF (debugging)**:
  * Adds missing gene/transcript IDs
  * Deduplicates features
  *	Fixes structural nesting issues in GTF entries (improper exon/transcript hierarchies), ensuring the overall annotation tree is coherent. 
  *	Splits overly large contigs when necessary (based on the genome statistics), accommodating tools that have contig size limits. The only mammalian species that I’ve seen need this is Monodelphis domestica.
*	**Shell script generation for references**: After cleaning, Cleanome can also write out shell scripts that will invoke the appropriate reference builder binary (e.g., the cellranger-arc CLI) with the cleaned GTF and accompanying FASTA, ready for indexing and full reference creation. 


### Installation

Make a anaconda environment with python>=3.6

```
conda install -c conda-forge ncbi-datasets-cli

conda install -c conda-forge -c bioconda ete3 gtfparse numpy pandas polars polars-lts-cpu pyarrow requests biopython tqdm 

git clone git@github.com:mtvector/cleanome.git
cd cleanome
pip install .
```

### Usage


To debug a gtf by adding missing gene and transcript fields, replacing missing gene fields with the gene_id, and other common issues in NCBI genome annotations:

```
debug_gtf file.gtf file.debug.gtf
```

See example_run.sh for an example script utilizing the full pipeline to download genomes, get statistics, debug gtfs and build cellranger-arc references.

Current pipeline functions include:
```
download_genomes --species_list ./species.txt --genome_dir ./genomes/
```
Download genomes from NCBI for all the species on the list. A utility for downloading ENSEMBL genomes for the list is included (have to write your own download_genomes.py for now)

```
get_genomes_and_stats --genome_dir ./genomes/ -o ./genome_info.csv -c
```
Collects the fastas and gtfs from all the genomes in a directory and calculates some simple statistics. 

```
make_cellranger_arc_sh --sh_scripts_dir ./submission_scripts/ --stats_csv ./genome_info.csv --output_dir ~/cellranger-arc --log_dir ~/log/ -cellranger_bin /path/to/cellranger-arc/bin/
```

Debugs gtfs (deduplicate transcripts/genes, add missing transcripts for exons and missing genes for transcripts, fills missing values with placeholders, split chromosomes that are too large, fix nesting)


