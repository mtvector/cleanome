# cleanome
Download genome annotations and prepare them to be made into reference annotations (e.g. with cellranger or cellranger-arc)


### Installation

Make a anaconda environment with python>=3.6

```
conda install -c conda-forge ncbi-datasets-cli```
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


