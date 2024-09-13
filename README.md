# DDSP (Define Secondary Structure of Proteins) Recreation Project

## About the Project

DPPP aims to simulate a very basic reproduction of the original DSSP software.

This README provides an overview of the functionality and usage of the dppp Python script.

### Getting Started

Prerequisites:
-The dppp_env.yaml file in the main repository.
        


### Installing

Use the following commands to setup your environment properly. Be sure to choose a suitable work directory first!

conda env create -f dppp_env.yaml

conda activate DPPP

### Usage

If not done in the previous step, use:
Activate environment: conda activate DPPP
to activate your environment.

Execute the main script, and make sure you have your PDB file on hand:

python script.py <PDB file path>

View Results:

    Custom Predictions: Results will be saved in custom_output.txt.
    DSSP Predictions: Results will be saved in real_output.txt.
    Memory Usage: Memory usage details will be printed in the console.
    Graphical Interface: A GUI window will display the secondary structure predictions.

You can reexecute this process with any PDB file. Once you're done, exit your environment with:

conda deactivate

Worry not, you can go back to it at any time!

### Performance Metrics

Time Measurement: The script measures the time taken for both the custom predictions and DSSP predictions.

Memory Usage: The script monitors memory usage and prints it upon completion.

### Built With

Python version 3.10.12
Libraries 
        -Native: Multiprocessing, Sys, Time, Tkinter
        -Additional: Biopython, 
### Contact Us
Project lead: wajih.bhh@outlook.com
Project link: https://github.com/WajihBHH/prog3

###Acknowledgements
https://pages.github.com/
