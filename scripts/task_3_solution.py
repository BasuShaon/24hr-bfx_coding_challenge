# %%
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations

# set up relative path
path = os.path.dirname('')
task_3_dir = os.path.join(path, '..', 'task_3')

# %% read in data

proteins = pd.read_csv(os.path.join(task_3_dir, 'proteins.txt'), names=['protein'])
compartments = pd.read_csv(os.path.join(task_3_dir, 'protein_compartments.csv'))  # 'protein_id', 'compartment_id'
interactions = pd.read_csv(os.path.join(task_3_dir, 'protein_interactions.txt'), sep=' ', header=None, names=['protein_A', 'protein_B'])

# rank order interactions
def reorder_columns(row):
    num_a = int(''.join(filter(str.isdigit, row['protein_A'])))
    num_b = int(''.join(filter(str.isdigit, row['protein_B'])))
    
    # Ensure the smaller value is in A and larger in B
    if num_a > num_b:
        row['protein_A'], row['protein_B'] = row['protein_B'], row['protein_A']
    return row

# Apply to each row
interactions_sorted = interactions.apply(reorder_columns, axis=1)

# Generate all possible pairs
all_pairs = pd.DataFrame(combinations(proteins['protein'], 2), columns=['protein_A', 'protein_B'])

#%%
#Identify all of the unique combinations of two proteins that do not reside in the same
#compartment, and also have not been previously identified to interact with each other

mask = all_pairs.set_index(['protein_A', 'protein_B']).index.isin(
    interactions_sorted.set_index(['protein_A', 'protein_B']).index
)

all_pairs_minus_interactions = all_pairs[~mask]
print(all_pairs_minus_interactions)

# %%
# Add compartment and component data to pairs
compartment_dict = compartments.set_index('protein_id')['compartment_id'].to_dict()
all_pairs['compartment_A'] = all_pairs['protein_A'].map(compartment_dict)
all_pairs['compartment_B'] = all_pairs['protein_B'].map(compartment_dict)

# Filter pairs: different compartment and different component (no direct or indirect interactions)
filtered_pairs = all_pairs[
    (all_pairs['compartment_A'] != all_pairs['compartment_B'])
]

print(filtered_pairs)

# %%
# Build connected components (manual approach, no NetworkX)
def find_components(interactions_df):
    parent = {}

    # find root component
    def find(u):
        if parent[u] != u:
            parent[u] = find(parent[u])  # Path compression
        return parent[u]

    # initialize each protein as its own parent
    all_proteins = set(interactions_df.protein_A).union(interactions_df.protein_B)
    for protein in all_proteins:
        parent[protein] = protein

    # union sets based on known interactions
    for _, row in interactions_df.iterrows():
        u, v = row['protein_A'], row['protein_B']
        u_root, v_root = find(u), find(v)
        if u_root != v_root:
            parent[v_root] = u_root

    # group proteins by their root parents
    components = {}
    for protein in all_proteins:
        root = find(protein)
        if root not in components:
            components[root] = set()
        components[root].add(protein)
    return components



# %%
# Get connected components

components = find_components(interactions)
# Map proteins to components
protein_component_map = {}
for component_id, proteins_set in enumerate(components.values()):
    for protein in proteins_set:
        protein_component_map[protein] = component_id

# %%

all_pairs['component_A'] = all_pairs['protein_A'].map(protein_component_map)
all_pairs['component_B'] = all_pairs['protein_B'].map(protein_component_map)

# Filter pairs: different compartment and different component (no direct or indirect interactions)
filtered_pairs = all_pairs[
    (all_pairs['compartment_A'] != all_pairs['compartment_B']) &
    (all_pairs['component_A'] != all_pairs['component_B'])
]