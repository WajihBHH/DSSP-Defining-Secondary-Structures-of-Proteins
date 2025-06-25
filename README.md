# DDSP (Define Secondary Structure of Proteins) Recreation Project

## Table of Contents
1. About the project and Installation
2. Getting PDB files
3. Running the code
4. Understanding the output
5. Analyzing multiple files
6. Built With
7. Contact us
8. Acknowledgements

## About the Project and Installation

This project aims to simulate a simpler reproduction of the original DSSP software.


It's focused on assigning alpha helix and beta sheet only.

### Install required packages


    pip install biopython numpy
    

### Optional : Install DSSP for reference comparison

    sudo apt-get install dssp # Ubuntu/Debian

    brew install dssp # macOS with Homebrew
    

## Getting PDB Files

### Direct download from RCSB PDB


#### Download 1ZAA (example protein)


    wget https://files.rcsb.org/download/1ZAA.pdb


#### Or download any other PDB file


    wget https://files.rcsb.org/download/[PDB_ID].pdb

### Manual download


Go to https://www.rcsb.org/


Search for your protein (e.g., "1ZAA")


Click "Download Files" → "PDB Format"


## Running the Code

### Through the command line

    python dssp_assignment.py 1zaa.pdb

### Or in the code

    pdb_file = "your_protein.pdb"
    
    Change this line then run python dssp_assignment.py

## Understanding the Output


### Output Files


[filename]_secondary_structure.txt : Tab-separated file with assignments


Contains columns: Residue, Chain, Number, Predicted


### Structure Type Codes


H : Alpha-helix

E : Beta-sheet

C : Coil/Random coil


### Metrics


Precision : Of predicted helices, how many are actually helices?


Recall : Of actual helices, how many did we predict correctly?


F1-Score : Balanced measure of precision and recall


### Hydrogen Bond Information


Shows donor → acceptor residue pairs


Energy in kcal/mol (more negative = stronger bond)


Distance in Angstroms

## Analyzing multiple files

--- Run this code

    import glob 
    
    from dssp_assignment import SecondaryStructureAssigner
    
    pdb_files = glob.glob("*.pdb")
    
    results = { }
    
    for pdb_file in pdb_files: 
    
        print(f"Analyzing { pdb_file } ...") 
        
        analyzer = SecondaryStructureAssigner(pdb_file) 
        
        assignments, accuracy, metrics = analyzer.run_analysis() 
        
        results[pdb_file] = accuracy 
        
        print ( f"Accuracy: {accuracy: .3f}\n")
        

## Built With

### Python version 3.10.12

### Libraries :
        math
        sys
        numpy
        Bio.PDB
        collections
        warnings

        
## Contact Us

Project lead: wajih.bhh@outlook.com

Project link: https://github.com/WajihBHH/prog3

## Acknowledgements

Kabch W., Sander C., et al. (1983)
