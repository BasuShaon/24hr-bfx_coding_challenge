# %%
import os
import pandas as pd
from itertools import combinations

class ProteinInteractionAnalyzer:
    def __init__(self, proteins_file, compartments_file, interactions_file):
        self.proteins = pd.read_csv(proteins_file, names=['protein'])
        self.compartments = pd.read_csv(compartments_file)
        self.interactions = pd.read_csv(interactions_file, sep=' ', header=None, names=['protein_A', 'protein_B'])
        self.interactions_sorted = self.rank_order_interactions()
        self.all_pairs = self.generate_all_pairs()
        self.components = self.find_connected_components()
        self.compartment_dict = self.compartments.set_index('protein_id')['compartment_id'].to_dict()
        self.protein_component_map = self.create_protein_component_map()

    def rank_order_interactions(self):
        interactions_sorted = self.interactions.copy()

        def reorder(row):
            num_a = int(''.join(filter(str.isdigit, row['protein_A'])))
            num_b = int(''.join(filter(str.isdigit, row['protein_B'])))
            if num_a > num_b:
                row['protein_A'], row['protein_B'] = row['protein_B'], row['protein_A']
            return row

        interactions_sorted = interactions_sorted.apply(reorder, axis=1)
        return interactions_sorted

    def generate_all_pairs(self):
        return pd.DataFrame(combinations(self.proteins['protein'], 2), columns=['protein_A', 'protein_B'])

    def filter_new_interactions(self):
        mask = self.all_pairs.set_index(['protein_A', 'protein_B']).index.isin(
            self.interactions_sorted.set_index(['protein_A', 'protein_B']).index
        )
        return self.all_pairs[~mask]

    def find_connected_components(self):
        parent = {}

        def find(u):
            if parent[u] != u:
                parent[u] = find(parent[u])
            return parent[u]

        all_proteins = set(self.interactions_sorted.protein_A).union(self.interactions_sorted.protein_B)
        for protein in all_proteins:
            parent[protein] = protein

        for _, row in self.interactions_sorted.iterrows():
            u, v = row['protein_A'], row['protein_B']
            u_root, v_root = find(u), find(v)
            if u_root != v_root:
                parent[v_root] = u_root

        components = {}
        for protein in all_proteins:
            root = find(protein)
            components.setdefault(root, set()).add(protein)

        return components

    def create_protein_component_map(self):
        protein_component_map = {}
        for component_id, proteins_set in enumerate(self.components.values()):
            for protein in proteins_set:
                protein_component_map[protein] = component_id
        return protein_component_map

    def filter_pairs_different_compartment_component(self):
        pairs = self.all_pairs.copy()
        pairs['compartment_A'] = pairs['protein_A'].map(self.compartment_dict)
        pairs['compartment_B'] = pairs['protein_B'].map(self.compartment_dict)
        pairs['component_A'] = pairs['protein_A'].map(self.protein_component_map)
        pairs['component_B'] = pairs['protein_B'].map(self.protein_component_map)

        return pairs[
            (pairs['compartment_A'] != pairs['compartment_B']) &
            (pairs['component_A'] != pairs['component_B'])
        ]


# Example usage:
path = os.path.dirname('')
task_3_dir = os.path.join(path, '..', 'task_3')

analyzer = ProteinInteractionAnalyzer(
    proteins_file=os.path.join(task_3_dir, 'proteins.txt'),
    compartments_file=os.path.join(task_3_dir, 'protein_compartments.csv'),
    interactions_file=os.path.join(task_3_dir, 'protein_interactions.txt')
)

new_interactions = analyzer.filter_new_interactions()
filtered_pairs = analyzer.filter_pairs_different_compartment_component()

print("New potential interactions:")
print(new_interactions)

print("\nFiltered pairs (different compartment and component):")
print(filtered_pairs)
