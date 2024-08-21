# cleanome
This is a python package that can be used to download genome annotations from multiple species and prepare them to be made into reference annotations (e.g. with cellranger or cellranger-arc). It is meant to standardize the process of debugging gtf files in order to make them compatible as gtfs provided by NCBI and ENSEMBL generally have one of a few problems that make them incompatible with creating genomics reference annotations (missing gene/transcript ids, duplicated genes or transcripts, contigs that are too large etc). The package has three functions: 1. Download fasta, gtf, and assembly metadata files for a list of species (by scientific name or taxid). 2. Generate statistics for the assemblies and debug them. 3. Write shell scripts for making [cellranger-arc] references.


### Installation

Make a anaconda environment with python>=3.6

```
conda install -c conda-forge ncbi-datasets-cli

conda install -c conda-forge -c bioconda ete3 gtfparse numpy pandas polars polars-lts-cpu pyarrow requests biopython
```

```
git clone git@github.com:mtvector/cleanome.git
cd cleanome
pip install .
```

### Usage

See example_run.sh for an example script to download genomes, get statistics, debug gtfs and build cellranger-arc references.

Current functions include:
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

Or to just debug a gtf:

```
debug_gtf file.gtf file.debug.gtf
```


