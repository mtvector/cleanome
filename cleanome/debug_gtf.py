#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Matthew Schmitz, Allen Institute, 2024

Sorts, deduplicates, fixes nesting, fills blanks, adds missing annotations for GTF files
usage:
python debug_gtf "{gtf_path}" "{gtf_output}"
"""


import pandas as pd
import gtfparse
import polars as pl
import os
import csv
import re
import sys
import tqdm
import numpy as np

def gtfparse_gtf_file(file_path):
    '''Read the GTF file into a Polars DataFrame'''
    return gtfparse.read_gtf(file_path)

def polars_to_pandas(df_polars):
    '''Convert Polars DataFrame to pandas DataFrame'''
    return df_polars.to_pandas()

def update_gene_coordinates_vectorized(df):
    # Only update rows with valid gene_id
    valid = df['gene_id'] != 'nan'
    # Ensure numeric types for the valid start and end values.
    df.loc[valid, 'start'] = pd.to_numeric(df.loc[valid, 'start'], errors='coerce')
    df.loc[valid, 'end'] = pd.to_numeric(df.loc[valid, 'end'], errors='coerce')
    # Compute min start and max end for each gene_id using a groupby transform.
    df['gene_min_start'] = df.groupby('gene_id')['start'].transform('min')
    df['gene_max_end'] = df.groupby('gene_id')['end'].transform('max')
    # Update gene rows: only modify rows where feature == 'gene' and gene_id is valid.
    mask = (df['feature'] == 'gene') & valid
    df.loc[mask, 'start'] = df.loc[mask, 'gene_min_start']
    df.loc[mask, 'end'] = df.loc[mask, 'gene_max_end']
    # Clean up helper columns.
    df.drop(columns=['gene_min_start', 'gene_max_end'], inplace=True)
    return df


def gtf_add_missing_features_optimized(df):
    def create_new_row(row, feature_type):
        new_row = row.copy()
        new_row['feature'] = feature_type
        new_row['frame'] = '.'
        new_row['score'] = '.'
        return new_row

    # Ensure IDs are strings for matching.
    df['transcript_id'] = df['transcript_id'].astype(str)
    df['gene_id'] = df['gene_id'].astype(str)

    existing_transcripts = set(df.loc[df['feature'] == 'transcript', 'transcript_id'])
    existing_genes = set(df.loc[df['feature'] == 'gene', 'gene_id'])

    new_transcript_rows = []
    new_gene_rows = []

    # Add missing transcript and gene rows from exon and transcript features.
    for feature in ['exon', 'transcript']:
        for _, row in tqdm.tqdm(df.loc[df['feature'] == feature].iterrows()):
            if feature == 'exon' and row['transcript_id'] not in existing_transcripts:
                new_transcript_rows.append(create_new_row(row, 'transcript'))
                existing_transcripts.add(row['transcript_id'])
            if row['gene_id'] not in existing_genes:
                new_gene_rows.append(create_new_row(row, 'gene'))
                existing_genes.add(row['gene_id'])

    new_df = pd.concat([df,
                        pd.DataFrame(new_transcript_rows),
                        pd.DataFrame(new_gene_rows)], ignore_index=True)
    new_df.loc[new_df['feature'] == 'exon', 'frame'] = '.'
    new_df.sort_values(by=['seqname', 'start'], inplace=True)

    # Use the vectorized update to adjust gene coordinates.
    new_df = update_gene_coordinates_vectorized(new_df)

    new_df = new_df.astype(str)
    new_df.fillna('nan', inplace=True)

    # Propagate gene/gene_name values within each gene_id group.
    for col in ['gene', 'gene_name']:
        if col in new_df.columns:
            new_df[col] = new_df[col].replace('nan', np.nan)
            new_df[col] = new_df[col].replace('', np.nan)
            new_df[col] = new_df.groupby('gene_id')[col].transform(lambda x: x.ffill().bfill())
            new_df[col] = new_df[col].fillna('nan')
    return new_df

def gtf_attribute_string_funfun(columns_to_condense):
    '''combines specified columns into gtf/gff2 attributes format'''
    def format_string(row):
        return '; '.join([f'{col} "{row[col]}"' for col in columns_to_condense if ((row[col]!='') and (row[col]!='nan'))])
    return format_string

def gtf_df_sort(df):
    '''sort a gtf df in the format required by cellranger, gene-transcript, exon nested and sorted'''
    def feature_sort(row):
        if row['feature'] == 'gene':
            return 1
        elif row['feature'] == 'transcript':
            return 2
        elif row['feature'] == 'exon':
            return 3
        else:
            return 4
    print('fix tx bug')
    df.loc[df['transcript_id']=='nan','transcript_id']=''
    df.loc[df['transcript_id'].isna(),'transcript_id']=''

    # Apply the custom sorting function
    df['feature_sort'] = df.apply(feature_sort, axis=1)
    
    # Sort by gene_id, feature_sort, and then start position
    df_sorted = df.sort_values(by=['gene_id','transcript_id','feature_sort','start','end'])
    
    # Drop the temporary sorting column
    df_sorted = df_sorted.drop('feature_sort', axis=1)
    return df_sorted

def write_gtf_df(df, output_file_path):
    '''write gtf df columns to gtf file format'''
    if len(df.columns)>8:
        if 'attribute' in df.columns:
            df.drop('attribute',axis=1,inplace=True)
        format_string=gtf_attribute_string_funfun(df.columns[8:])
        df.loc[:,'attribute']=df.apply(format_string, axis=1)
    columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    #####FOR SOME REASON GTF MUST BE IN ASCII ENCODING FOR CELLRANGER!!!!!!!
    df.loc[:,columns].to_csv(output_file_path, sep='\t', index=False, header=False,quoting=csv.QUOTE_NONE,encoding='ascii')

def make_unique(df, column, feature_name='gene'):
    """
    Append a suffix to duplicate values in the specified column of a DataFrame
    to make them unique. Only applies to rows where feature == "gene".
    """
    if 'feature' not in df or column not in df:
        return df  # Safety check

    # Filter to process only rows where feature is 'gene'
    genes = df[df['feature'] == feature_name]

    # Find duplicates and make them unique
    duplicates = genes.duplicated(subset=[column], keep=False)
    unique_counter = {}  # To keep track of how many times we've seen a value
    for i, is_duplicate in enumerate(duplicates):
        if is_duplicate:
            value = genes.iloc[i][column]
            count = unique_counter.get(value, 0)
            # Update the DataFrame with the new unique value
            df.loc[genes.index[i], column] = f"{value}__{count+1}"
            unique_counter[value] = count + 1

    return df

def deduplicate_gtf(pandas_df,extra_allowed_cols=[]):
    allowed_cols=['seqname','source','feature','start','end','score','strand','frame','gene_id','transcript_id','exon_id','biotype', 'db_xref', 'description', 'gbkey','gene','gene_symbol','hgnc_symbol','gene_name']+extra_allowed_cols
    pandas_df=pandas_df.loc[:,pandas_df.columns.isin(allowed_cols)]
    pandas_df=pandas_df.astype(str).fillna('nan')
    pandas_df=pandas_df.replace(r'^\s*$', 'nan', regex=True)
    pandas_df=make_unique(pandas_df,'gene_id')
    pandas_df=make_unique(pandas_df,'gene')
    pandas_df=make_unique(pandas_df,'transcript_id',feature_name='transcript')
    return pandas_df

def fill_missing_names_with_id(pandas_df,name_column='gene',fill_column='gene_id'):
    missing_names = (pandas_df[name_column] == 'nan') | (pandas_df[name_column] == '') | (pandas_df[name_column].isna())
    pandas_df.loc[missing_names,name_column] = pandas_df.loc[missing_names,fill_column]
    return pandas_df

def transfer_gtf_header(source_gtf, target_gtf, output_gtf):
    header_lines = []

    with open(source_gtf, 'r') as src_file:
        for line in src_file:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                break

    with open(target_gtf, 'r') as tgt_file:
        target_content = tgt_file.readlines()

    # Write the header lines and target content to a new output file
    with open(output_gtf, 'w') as out_file:
        out_file.writelines(header_lines)
        out_file.writelines(target_content)

def main():
    file_path=sys.argv[1]
    polars_df = gtfparse_gtf_file(file_path)
    pandas_df = polars_to_pandas(polars_df)
    pandas_df = gtf_df_sort(pandas_df)
    df_with_transcripts = gtf_add_missing_features_optimized(pandas_df)
    df_with_transcripts = fill_missing_names_with_id(df_with_transcripts) if len(sys.argv)<4 else fill_missing_names_with_id(df_with_transcripts,sys.argv[3])#,sys.argv[3]
    print('added features')
    df_with_transcripts=df_with_transcripts.loc[~df_with_transcripts['transcript_id'].astype(str).str.contains('unknown'),:]
    df_with_transcripts=deduplicate_gtf(df_with_transcripts)
    df_with_transcripts = gtf_df_sort(df_with_transcripts)
    print('sorted')
    df_with_transcripts['score']='.'#the score column is worthless for cellranger so who cares
    if 'attribute' in df_with_transcripts.columns:
        df_with_transcripts.drop('attribute',axis=1,inplace=True)
    write_gtf_df(df_with_transcripts, sys.argv[2])
    transfer_gtf_header(sys.argv[1], sys.argv[2], sys.argv[2])
    
if __name__ == "__main__":
    main()
