# Relation coding challenge (Shaon)

This repository contains solutions for the 3 tasks. 
The project follows the following repository structure:

## Project Structure

```
project_root/           
│-- task_1/             
│   ├── human_proteins_dirty.csv
│   ├── human_proteins_clean.csv # cleaned file from code
│-- task_2/             
│   ├── approved_drugs.csv
│   ├── spike_positive_compounds.csv
│   ├── stat_report.csv # feature-wise stat report from code
│-- task_3/         
│   ├── proteins.csv
│   ├── protein_compartments.csv
│   ├── protein_interactions.csv
│
│-- scripts/             
│   │-- task_1_solution.py
│   │-- task_2_solution.pynb
│   │-- task_3_solution.py
│   │-- utils.py        
│
│-- requirements.txt     
│-- README.md           
```

## Setup Instructions

### 1. Install Dependencies
Ensure you have Python installed, then run:
```sh
cd Relation_coding_challenge
pip install -r requirements.txt
```

### 2. Running the Scripts for Task 1 and 3
Navigate to the `scripts` folder and run the respective script for each task:
```sh
python scripts/task_1_solution.py
python scripts/task_3_solution.py
```

### 3. Running the Jupyter Notebook for Task 2
Launch Jupyter Notebook:
```sh
jupyter notebook
```
Once the notebook interface opens, navigate to the notebook/ folder and open task_1_solution.ipynb.


