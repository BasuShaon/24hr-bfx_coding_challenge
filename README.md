# Relation coding challenge (Shaon)

This repository contains solutions for the 3 tasks. 
The project follows the following repository structure:

## Project Structure

```
project_root/           
│-- task_1/             
│   ├── human_proteins_dirty.csv
│   ├── human_proteins_clean.csv    # cleaned file for task_1
│   ├── human_proteins_clean_patternX.csv     

│-- task_2/             
│   ├── approved_drugs.csv
│   ├── spike_positive_compounds.csv
│   ├── stat_report.csv             # stat report for task_2
│
│-- task_3/         
│   ├── proteins.csv
│   ├── protein_compartments.csv
│   ├── protein_interactions.csv
│   ├── protein_answer1.csv         # PPIs from task_3 (Q1 criteria)
│   ├── protein_answer2.csv         # PPIs from task_3 (Q2 criteria)
│
│-- scripts/             
│   │-- task_1_solution.py
│   │-- task_2_solution.ipynb
│   │-- task_3_solution.py
│   │-- protein_network_analyzer.py        
│
│-- requirements.txt     
│-- README.md           
```

## Setup Instructions

### 1. Install Dependencies
Ensure you have **Python 3.11+** installed. Then, navigate to the root of the repository and install dependencies:
```sh
cd Relation_coding_challenge
pip install -r requirements.txt
```

### 2. Running Scripts for Task 1 and 3
Navigate to the `scripts` folder and run the respective script for each task:
```sh
Python3 scripts/task_1_solution.py
Python3 scripts/task_3_solution.py
```

### 3. Running Jupyter Notebook for Task 2
Launch Jupyter Notebook:
```sh
jupyter notebook
```
Once the notebook interface opens, navigate to the `scripts` folder and open `task_1_solution.ipynb`.


