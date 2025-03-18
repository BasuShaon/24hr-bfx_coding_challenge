# %%
import os
import numpy as np
import pandas as pd
from itertools import product

# Get relative path for data directory
current_path = os.path.dirname(__file__)
data_dir = os.path.join(current_path, '..', 'task_1')

def diagnose_corruptions(file_path):
    """
    Reads a potentially corrupt CSV file and identifies warnings/issues using pandas.read_csv()
    on_bad_lines = 'warn' functionality 

    Parameters
    ----------
    data : str
        relative filepath to corrupt data

    """
    print('\nDiagnose .csv corruptions by index (rows):')
    print(pd.read_csv(file_path, on_bad_lines='warn', index_col=0))


def clean_protein_file(input_path, output_path):
    """
    Cleans the corrupted human protein file by fixing missing values and saves a corrected version.
   
    Parameters
    ----------
    input_path : str
        relative filepath to corrupt data (ending with .csv)
    output_path : str
        relative filepath to save cleaned data (ending with .csv)
    """
    # Read the corrupt file as 4 fields to acount for the 3 - 4 field mismatch, skipping header row
    protein_df = pd.read_csv(input_path, names=[0, 1, 2, 3], skiprows=1, index_col=0)

    # Identify rows where column 2 has missing values (corruption)
    corrupt_rows = protein_df.loc[protein_df[2].isna()].index

    # Fix the corropted rows: move values from column 3 to column 2
    protein_df.loc[corrupt_rows, 1] = protein_df.loc[corrupt_rows, 2].values

    # Save cleaned data, keeping only relevant columns (the expected 3 fields)
    protein_df.iloc[:, 0:2].to_csv(output_path, index=True, header=['protein_id', 'sequence'])

    print(f"\nCleaned file saved to: {output_path}")

def find_pattern_combinations(df): 
    """
    Method 1 to search for matching possible hits (uses combinatorics):
    Uses combinatorial approach to generate sequences matching T-AV-T-*-T motif / pattern, 
    then performs search for protein sequences in the dataset that contain the possible patterns

    Parameters
    ----------
    df : pd.DataFrame
        The cleaned protein_id, sequence dataframe

    Returns
    -------
    matches : pd.DataFrame
        dataframe containing protein hits with matches to pattern X (T-AV-T-*-T motif)

    """
    # Possible values for the ordered amino acid entries in pattern X (the query motif) 
    entry_1 = ['T']
    entry_2 = ['A', 'V']
    entry_3 = ['T']
    entry_4 = list('ACDEFGHIKLMNPQRSTVWY')  
    entry_5 = ['T']

    # Generate all combinations that match the pattern (1 x 2 x 1 x 20 x 1) = 40 
    combinations = [''.join(p) for p in product(entry_1, entry_2, entry_3, entry_4, entry_5)]
    print(f'\n{len(combinations)} possible motifs with pattern X combination:', combinations)
    # Convert list into regex pattern (replace white space with '|')
    pattern = '|'.join(combinations)  

    # Find matching motifs in proteins in the dataframe
    matches = df[df.iloc[:, 1].str.contains(pattern, na=False)]
    matches.set_index('protein_id', inplace=True)

    print('\nIdentified hits using combinatorics:\n', matches)
    return matches

def find_pattern_regex(df):
    """
    Method 2 to search for matching possible hits (uses regex):
    Uses regex to identify protein sequences in the dataset that contain
    T-AV-T-*-T motif / pattern. 

    Parameters
    ----------
    df : pd.DataFrame
        The cleaned protein_id, sequence dataframe

    Returns
    -------
    matches : pd.DataFrame
        dataframe containing protein hits with matches to pattern X (T-AV-T-*-T motif)

    """
    regex_pattern = r'T[AV]T[A-Z]T'  
    
    matches = df[df.iloc[:, 1].str.contains(regex_pattern, regex=True, na=False)]
    matches.set_index('protein_id', inplace=True)

    print('\nIdentified hits using regex:\n', matches)
    return matches

# Set filepaths
dirty_file = os.path.join(data_dir, 'human_proteins_dirty.csv')
clean_file = os.path.join(data_dir, 'human_proteins_clean.csv')

# Step 1: Diagnose corruptions
corrupt_df = diagnose_corruptions(dirty_file)

# Step 2: Clean the file and save a corrected version
clean_protein_file(dirty_file, clean_file)

# Step 3: Load the cleaned file for analysis
protein_df = pd.read_csv(clean_file, index_col=0)

# Step 4: Find sequences matching the T-AV-T-*-T motif pattern
matches_combinatorics = find_pattern_combinations(protein_df)
matches_regex = find_pattern_regex(protein_df)

# Step 5: Save output
matches_combinatorics.to_csv(os.path.join(data_dir, 'human_proteins_clean_patternX.csv'))
print(f"\nHits saved to: {os.path.join(data_dir, 'human_proteins_clean_patternX.csv')}")
