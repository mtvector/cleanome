"""
Author: Matthew Schmitz, Allen Institute, 2024
"""
import pandas as pd
import numpy as np
from Bio imort SeqIO
import csv

def zcumsum(iterable):
    cumulative_sum=[0]
    total=0
    for item in iterable:
        total += item
        cumulative_sum.append(total)
    cumulative_sum=cumulative_sum+[float("inf")]
    return cumulative_sum


class ChromosomeSplitter:
    """
    Splits contigs greater than a threshold into multiple contigs, and corrects associated GTF file

    Detecting breakpoints is as follows: 
    for contigs larger than the threshold, identify the theoretical 
    split indices that will split it into n equal size contigs 
    smaller than the threshold. Then for each split point, 
    look at the nearest 50 intergenic regions (using the gtf gene annotations) 
    and pick the midpoint of the largest intergenic region as the 
    optimal split point for each split.

    Takes a fasta and gtf and Can also take a dictionary of lists of breakpoints as input to skip detection step.
    ##detecting new from scratch
    ChromosomeSplitter(fasta_file,gtf_file,new_fasta_file,new_gtf_file)
    
    # Example usage of using specified breakpoints
    import os
    top_path='/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/genomes/Opossum/raw_files/'
    fasta_file=os.path.join(top_path,'GCF_027887165.1_mMonDom1.pri_genomic.fna')
    new_fasta_file=os.path.join(top_path,'GCF_027887165.1_mMonDom1.pri_genomic_split.fna')
    gtf_file=os.path.join(top_path,'GCF_027887165.1_mMonDom1.pri_genomic.gtf')
    new_gtf_file=os.path.join(top_path,'GCF_027887165.1_mMonDom1.pri_genomic_split.gtf')
    
    split_locations={
        'NC_077227.1': [495286654],
        'NC_077228.1': [502619815]
    }
    
    new_sequences=ChromosomeSplitter.read_fasta_and_split(fasta_file, split_locations)
    ChromosomeSplitter.write_new_fasta(new_sequences, new_fasta_file)
    ChromosomeSplitter.gtf_chrom_name_change(gtf_file, split_locations, new_gtf_file)
    """
    def __init__(self,fasta_file,gtf_file,new_fasta_file,new_gtf_file,length_threshold=5e8,split_locations=None):
        if split_locations is None:
            split_locations=self.find_split_points(fasta_file, gtf_file, length_threshold)

        self.split_locations=split_locations
        new_sequences=self.read_fasta_and_split(fasta_file, split_locations)
        self.write_new_fasta(new_sequences, new_fasta_file)
        self.gtf_chrom_name_change(gtf_file, split_locations, new_gtf_file)


    @staticmethod
    def read_fasta(fasta_file):
        sequences = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences[record.id] = len(record.seq)
        return sequences

    @staticmethod
    def read_gtf_as_bed(gtf_file,gene_only=True):
        # Read the GTF file
        df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None,
                         names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    
        # Adjust start position for BED format (0-based)
        df['start'] = df['start'] - 1
    
        # Extract gene name or other identifiers from the attributes column
        df['name'] = df['attribute'].str.extract('gene_id "([^"]+)"')
    
        # If gene name is not found, extract gene_id
        df.loc[df['name'].isna(),'name']=df['attribute'].str.extract('gene_name "([^"]+)"')[df['name'].isna()].to_numpy()
       
        df['extra'] = df['attribute']
        if gene_only:
            df=df.loc[df['feature']=='gene',:]
            
        bed_df = df[['seqname', 'start', 'end', 'name', 'score', 'strand', 'extra']]   
        return bed_df

    @staticmethod
    def find_nearest_intergenic_regions(genes, split_point, num_regions=50):
        # Sort genes based on the distance to the split point
        genes.loc[:,'distance'] = genes['start'].subtract(split_point).abs()
        nearest_genes = genes.sort_values(by='distance').head(num_regions)
        nearest_genes = nearest_genes.sort_values(by='start')
        all_sites=np.array(list(nearest_genes['start'].to_numpy())+list(nearest_genes['end'].to_numpy()))
        # Find intergenic regions among these genes
        intergenic_regions = []
        prev_end = nearest_genes['end'].to_numpy()[0]
        for _, row in nearest_genes.iterrows():
            #First is normal intergenic, second checks against possibility of nested genes
            if (row['start'] > prev_end) and (row['start']==all_sites[all_sites>prev_end].min()):
                intergenic_regions.append((prev_end, row['start']))
            prev_end = row['end']
        return intergenic_regions

    @classmethod
    def find_split_points(cls,fasta_file, gtf_file, length_threshold):
        sequences = cls.read_fasta(fasta_file)
        genes = cls.read_gtf_as_bed(gtf_file)
    
        split_points = {}
    
        for contig, length in sequences.items():
            if length <= length_threshold:
                continue
    
            # Calculate the number of segments needed
            num_segments = int(np.ceil(length / length_threshold))
            # Identify theoretical split indices
            theoretical_splits = [i * length // num_segments for i in range(1, num_segments)]
            # Filter genes on this contig
            contig_genes = genes[genes['seqname'] == contig]
            optimal_splits = []
            for split in theoretical_splits:
                # Find nearest intergenic regions to this split
                intergenic_regions = cls.find_nearest_intergenic_regions(contig_genes, split)
                # Choose the midpoint of the largest intergenic region
                largest_intergenic = max(intergenic_regions, key=lambda x: x[1] - x[0])
                optimal_split = (largest_intergenic[0] + largest_intergenic[1]) // 2
                optimal_splits.append(optimal_split)
    
            split_points[contig] = optimal_splits
    
        return split_points
    
    @staticmethod
    def read_fasta_and_split(file_path, split_locations):
        """
        Reads in a fasta and split chromosomes at points specicified indicated in split_locations dict
        """
        new_sequences={}
        for record in SeqIO.parse(file_path, "fasta"):
            chrom=record.id
            if chrom in split_locations:
                starts=[0] + split_locations[chrom]
                ends=split_locations[chrom] + [len(record.seq)]
                for i, (start, end) in enumerate(zip(starts, ends)):
                    new_id=f"{chrom}__part{i+1}"
                    new_sequences[new_id]=record.seq[start:end]
            else:
                new_sequences[chrom]=record.seq
        return new_sequences

    @staticmethod
    def write_new_fasta(new_sequences, output_file):
        """
        Writes a seqio genome object to fasta
        """
        with open(output_file, "w") as output_handle:
            for chrom, seq in new_sequences.items():
                SeqIO.write(SeqIO.SeqRecord(seq, id=chrom, description=""), output_handle, "fasta")
    
    @staticmethod
    def gtf_chrom_name_change(gtf_file, split_locations, output_file):
        """
        Replaces split chromosome names in associated gtf and corrects their start and end point
        """
        with open(gtf_file) as gtf_handle, open(output_file, "w") as out_handle:
            gtf_reader=csv.reader(gtf_handle, delimiter="\t")
            gtf_writer=csv.writer(out_handle, delimiter="\t")
            split_starts={x:zcumsum(split_locations[x]) for x in split_locations.keys()}
            
            for row in gtf_reader:
                if row[0].startswith('#'):
                    gtf_writer.writerow(row)
                    continue
                chrom=row[0]
                start=int(row[3])
                end=int(row[4])
    
                if chrom in split_locations.keys():
                    for idx in range(len(split_starts[chrom])-1):
                        if start < split_starts[chrom][idx+1] and start > split_starts[chrom][idx]:
                            new_chrom_name=f'{chrom}__part{idx+1}'
                            row[0]=new_chrom_name
                            row[3]=str(start - split_starts[chrom][idx])
                            row[4]=str(end - split_starts[chrom][idx])
                            break
    
                gtf_writer.writerow(row)
