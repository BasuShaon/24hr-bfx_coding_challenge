# %%
import os
import sys
import importlib
sys.path.append(os.path.join(os.path.dirname(__file__)))
import protein_network_analyzer 

# %%
def run_analysis():

    # Import classes and reload network analyzer, set directory path for task 3
    importlib.reload(protein_network_analyzer)
    path = os.path.dirname(__file__)
    task_3_dir = os.path.join(path, '..', 'task_3')

    # Initialize network analyzer class
    analyzer = protein_network_analyzer.ProteinNetworkAnalyzer(
        proteins_file=os.path.join(task_3_dir, 'proteins.txt'),
        compartments_file=os.path.join(task_3_dir, 'protein_compartments.csv'),
        interactions_file=os.path.join(task_3_dir, 'protein_interactions.txt')
    )

    # Answer question 1
    new_interactions = analyzer.select_crosscompartment_unobserved_interactions()
    print("New potential interactions (unobserved and different compartments):")
    print(new_interactions)
    new_interactions.to_csv(os.path.join(task_3_dir, 'protein_question1.csv'))
    print(f"\nQ1 answer saved to: {os.path.join(task_3_dir, 'proteins_question1.csv')}")

    # Answer question 2
    filtered_pairs = analyzer.select_crossnetwork_crosscompartment_interactions()
    print("\nFiltered pairs (different compartment and different network):")
    print(filtered_pairs)
    filtered_pairs.to_csv(os.path.join(task_3_dir, 'protein_question2.csv'))
    print(f"\nQ2 answer saved to: {os.path.join(task_3_dir, 'proteins_question2.csv')}")


if __name__ == "__main__":
    run_analysis()

