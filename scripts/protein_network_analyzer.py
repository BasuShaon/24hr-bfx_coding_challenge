# %%
import os
import pandas as pd
from itertools import combinations

class ProteinNetworkAnalyzer:
    """
    A class to analyze protein-protein interactions by processing protein data,
    compartments, and known interactions. It identifies new in-direct interactions
    and filters them based on biological constraints such as compartments and
    connectivity in interaction networks.

    Attributes
    ----------
    proteins : pd.DataFrame
        List of proteins

    compartments : pd.DataFrame 
        Mapping of proteins to cell compartments
    
    compartment_dict : dict
        Maps protein IDs to their cell compartments

    interactions : pd.DataFrame
        Known protein interactions  

    interactions_sorted : pd.DataFrame
        Interactions sorted in a standard order

    all_pairs : pd.DataFrame
        All possible unique protein pairs

    networks : dict
        Connected networks of the interaction network

    protein_network_map : dict
        Maps proteins to their connected networks

    """

    def __init__(self, proteins_file, compartments_file, interactions_file):
        """
        Initializes the class 
        
        Parameters
        ----------
        proteins_file : str
            Path to the protein id file

        compartments_file : str
            Path to the cell compartments metadata file
            
        interactions_file : str
            Path to the pair-wise direct interactions file

        """

        # load all proteins into a DataFrame
        self.proteins = pd.read_csv(proteins_file, names=['protein'])
        # load proteins and there annotated cell compartment into a DataFrame
        self.compartments = pd.read_csv(compartments_file)
        self.compartment_dict = self.compartments.set_index('protein_id')['compartment_id'].to_dict()
        # load the observed protein-protein interactions into a DataFrame
        self.interactions = pd.read_csv(interactions_file, sep=' ', header=None, names=['protein_A', 'protein_B'])
        # sort the observed protein-protein interactions so they are ordered by increasing protein ID (numerical)
        self.interactions_sorted = self.rank_order_interactions()
        # generate all theoretically possible protein-protein interactions
        self.all_pairs = self.generate_all_pairs()
        # find indirect-interaction network on observed protein-protein interactions using Union-Find algorithm
        self.networks = self.find_connected_networks()
        # map indirect-interactions into subgraph networks
        self.protein_network_map = self.create_protein_network_map()

    def rank_order_interactions(self):
        """
        Sorts interactions so that protein_A always has the smaller numerical ID.

        Returns
        ----------
        interactions_sorted : pd.DataFrame
            DataFrame of sorted interactions.

        """
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
        """
        Generates all possible unique protein-protein pairs.

        Returns
        ----------
        all_pairs : pd.DataFrame
            DataFrame containing all possible pairs.
            
        """
        all_pairs = pd.DataFrame(combinations(self.proteins['protein'], 2), columns=['protein_A', 'protein_B'])
        return all_pairs

    def select_crosscompartment_unobserved_interactions(self):
        """
        Identifies interaction pairs that are cross-cell-ompartment and 
        outside the known (unobserved) interactions set

        Returns
        ----------
        all_pairs_cc_no_observed : pd.DataFrame
            DataFrame containing only unobserved interactions.

        """
        pairs = self.all_pairs.copy()
        pairs['compartment_A'] = pairs['protein_A'].map(self.compartment_dict)
        pairs['compartment_B'] = pairs['protein_B'].map(self.compartment_dict)

        # cross-compartment interactions selection
        all_pairs_cc = pairs[
            (pairs['compartment_A'] != pairs['compartment_B'])]
        
        # mask for observed interactions
        mask = all_pairs_cc.set_index(['protein_A', 'protein_B']).index.isin(
            self.interactions_sorted.set_index(['protein_A', 'protein_B']).index)

        # unobserved interactions selection (inverse mask indexing on cross-compartmartment interactions df)
        all_pairs_cc_no_observed = all_pairs_cc[~mask]

        return all_pairs_cc_no_observed

    def find_connected_networks(self):
        """
        Main method which identifies connected network in the protein interaction network using the Union-Find algorithm.

        This method constructs a disjoint-set (Union-Find) structure to determine which proteins are
        directly or indirectly connected based on known interactions. Each connected network represents
        a group of proteins that interact either directly or through a chain of interactions. Each network
        forms a non-joint set (i.e. connection in venn diagram viz), and thus not connected by edges 
        (protein-protein interaction).

        Returns
        ----------
        networks : dict
            A dictionary where,
            - Keys are the main root protein (parent) of a connected network.
            - Values are sets of proteins belonging to that network.

        """
        main_root = {}

        def find(u):
            """
            Finds the main root (parent) of the connected network that protein `u` belongs to.
            
            This function implements **path compression** to optimize Union-Find operations.
            
            Parameters
            ----------
            u : str 
                The protein ID
            
            Returns
            ----------
            main_root : str
                The root protein ID representing the connected network
                
            """
            if main_root[u] != u:
                main_root[u] = find(main_root[u])
            return main_root[u]

        # initialize all proteins as their own roots (self-loops initially)
        all_proteins = set(self.interactions_sorted.protein_A).union(self.interactions_sorted.protein_B)
        for protein in all_proteins:
            main_root[protein] = protein  

        # merge network using Union-Find
        for _, row in self.interactions_sorted.iterrows():
            u, v = row['protein_A'], row['protein_B']
            u_root, v_root = find(u), find(v)
            if u_root != v_root:
                main_root[v_root] = u_root

        # group proteins by their representative network root
        networks = {}
        for protein in all_proteins:
            root = find(protein)
            networks.setdefault(root, set()).add(protein)

        return networks

    def create_protein_network_map(self):
        """
        Maps each protein to its corresponding connected network ID.

        Returns
        ----------
        protein_network_map : dict
            A dictionary mapping proteins to network IDs.

        """

        protein_network_map = {}
        for network_id, proteins_set in enumerate(self.networks.values()):
            for protein in proteins_set:
                protein_network_map[protein] = network_id
        return protein_network_map

    def select_crossnetwork_crosscompartment_interactions(self):
        """
        Filters protein pairs that belong to different compartments and different indirect-networks.
        Therefore, only selects for proteins that should not interact (different networks and interaction networks)

        Returns
        ----------
        pairs : pd.DataFrame
            DataFrame containing filtered protein pairs.
            
        """
        
        pairs = self.all_pairs.copy()
        pairs['compartment_A'] = pairs['protein_A'].map(self.compartment_dict)
        pairs['compartment_B'] = pairs['protein_B'].map(self.compartment_dict)
        pairs['network_A'] = pairs['protein_A'].map(self.protein_network_map)
        pairs['network_B'] = pairs['protein_B'].map(self.protein_network_map)

        return pairs[
            (pairs['network_A'] != pairs['network_B']) &
            (pairs['compartment_A'] != pairs['compartment_B']) 
        ]

