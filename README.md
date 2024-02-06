# cleanome
Download genome annotations and prepare them to be made into reference annotations (e.g. with cellranger or cellranger-arc)

### Usage
Current functions include:
```
cleanome/download_genomes.py --species_list ./species.txt --genome_dir ./genomes/ --stats_csv ./stats.csv
```
Download genomes from NCBI for all the species on the list. A utility for downloading ENSEMBL genomes for the list is included (have to write your own download_genomes.py for now)

```
cleanome/make_cellranger_arc_sh.py --sh_scripts_dir ./submission_scripts/ ---debug_script cleanome/debug_gtf.py --output_dir ~/cellranger-arc --log_dir ~/log/ -cellranger_bin /path/to/cellranger-arc/bin/
```
Debugs gtfs (deduplicate transcripts/genes, add missing transcripts for exons and missing genes for transcripts, fills missing values with placeholders, split chromosomes that are too large, fix nesting)

