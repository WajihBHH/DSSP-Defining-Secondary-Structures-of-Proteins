import sys
import time
import tkinter as tk

from Bio.PDB import PDBParser, PPBuilder, DSSP
from tkinter import ttk
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from memory_profiler import memory_usage

class Dppp:
    """
    A class to emulate the historic DSSP method as a secondary structure predictor.
    """

    def __init__(self, file_path):
        self.file_path = file_path


    def read_pdb_file(self):
        """Read PDB file and return structure."""

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein_structure', self.file_path)

        return structure


    def get_structure_info(self, structure):
        """Print information about the PDB structure."""

        num_models = len(structure)
        num_chains = num_residues = num_atoms = 0

        for model in structure:
            for chain in model:
                num_chains += 1
                for residue in chain:
                    num_residues += 1
                    for atom in residue:
                        num_atoms += 1

        print(f"Number of models: {num_models}")
        print(f"Number of chains: {num_chains}")
        print(f"Number of residues: {num_residues}")
        print(f"Number of atoms: {num_atoms}")


    def predict_multiprocessing(self, structure):
        """
        Predict secondary structure elements (a-helix and b-strands) using multiprocessing.
        """

        # Measure time needed.
        start_time = time.time()

        residues_data = []

        for model in structure:
            for chain in model:
                residues = list(chain)
                for i, residue_i in enumerate(residues):
                    residues_data.append((residue_i, i, residues))

        # Start the multiprocessing.
        with Pool(processes=cpu_count()) as pool:
            # Parallelize prediction
            results = list(tqdm(pool.imap(self._predict_for_residue, residues_data), \
            total=len(residues_data), desc="Squeezing PC juice...", unit="residu"))

        # Send predictions into an output file.
        predictions = [(i, prediction) for i, prediction in results]
        with open("custom_output.txt", "w") as out:
            print(f"Prediction Table \n", predictions, file=out)

        # End time measure and print result.
        end_time = time.time()
        dpppTime = end_time - start_time
        print(f"Time taken for DPPP prediction: {dpppTime:.2f} seconds.")

        return predictions


    def _predict_for_residue(self, data):
        residue_i, i, residues = data
        predictions = []
        prediction = '-'

        # Calculate a-helix by checking hydrogen bond between i and i+4.
        if i + 4 < len(residues):
            for j in range(i + 4, min(i + 10, len(residues))):
                residue_j = residues[j]
                energy = self.calculate_hydrogen_bond_energy(residue_i, residue_j)
                if energy < -0.5:  # Adjusted threshold for more sensitivity
                    prediction = 'H'
                    break

        # Calculate b-strand interactions. (parallel/antiparallel b-strands)
        if prediction == '-':  # Only check for b-strands if not already predicted as alpha helix
            for j in range(i + 1, len(residues)):
                residue_j = residues[j]
                energy = self.calculate_hydrogen_bond_energy(residue_i, residue_j)
                if energy < -0.6 and abs(i - j) > 3:
                    prediction = 'E'
                    break

        # Check for pi helices 
        if prediction == '-' and i + 4 < len(residues):
            for j in range(i + 4, min(i + 8, len(residues))):
                residue_j = residues[j]
                energy = self.calculate_hydrogen_bond_energy(residue_i, residue_j)
                if energy < -0.4:
                    prediction = 'I'
                    break

        # Check for 3/10 helices
        if prediction == '-' and i + 3 < len(residues):
            for j in range(i + 3, min(i + 6, len(residues))):
                residue_j = residues[j]
                energy = self.calculate_hydrogen_bond_energy(residue_i, residue_j)
                if energy < -0.3:
                    prediction = 'G'
                    break

        return i, prediction


    def calculate_hydrogen_bond_energy(self, residue1, residue2):
        """
        Calculate the hydrogen bond energy between two residues.
        """

        # Check if necessary atoms exist in both residues.
        if residue1.get_resname() not in ['ARG', 'LYS', 'HIS', 'SER', 'THR'] \
        or \
        residue2.get_resname() not in ['ASP', 'GLU']:
            return float('inf')

        # Get atom coordinates (O, N, C) from the residues.
        O1 = residue1['O'].get_vector() if 'O' in residue1 else None
        N2 = residue2['N'].get_vector() if 'N' in residue2 else None
        C1 = residue1['C'].get_vector() if 'C' in residue1 else None

        if O1 is None or N2 is None or C1 is None:
            return float('inf')

        # Calculate the distances between atoms.
        r_ON = (O1 - N2).norm()
        r_CN = (C1 - N2).norm()
        r_OC = (O1 - C1).norm()

        # Calculate hydrogen bond energy using a simplified formula.
        energy = 0.084 * (1 / r_ON + 1 / r_CN - 1 / r_OC) * 332

        # Print energy values, mostly for debugging purposes.
        #print(f"Residue pair: {residue1.get_resname()}-{residue2.get_resname()}")
        #print(f"Distances: r_ON={r_ON}, r_CN={r_CN}, r_OC={r_OC}")
        #print(f"Energy: {energy}")
        return energy


    # Calculation including Hydrogen. Less accurate and optional.

    #     def calculate_hydrogen_bond_energy(self, residue1, residue2):
    #         """
    #         Calculate the hydrogen bond energy between two residues.
    #         """
    #         if 'CA' not in residue1 or 'N' not in residue2 or 'C' not in residue1 or 'H' not in residue2:
    #             return float('inf')

    #         try:
    #             O1 = residue1['O'].get_vector()
    #             N2 = residue2['N'].get_vector()
    #             C1 = residue1['C'].get_vector()
    #             H2 = residue2['H'].get_vector()
    #         except KeyError:
    #             return float('inf')

    #         r_ON = (O1 - N2).norm()
    #         r_CH = (C1 - H2).norm()
    #         r_OH = (O1 - H2).norm()
    #         r_CN = (C1 - N2).norm()

    #         energy = 0.084 * (1 / r_ON + 1 / r_CH - 1 / r_OH - 1 / r_CN) * 332
    #         return energy


    def calculate_dssp_predictions(self, structure):
        """
        Use DSSP to get secondary structure predictions for comparison.
        """

        # Measure time needed.
        start_time = time.time()

        dssp_predictions = []
        for model in structure:
            for chain in model:
                dssp_obj = DSSP(model, self.file_path)
                for res in dssp_obj:
                    res_id, dssp_code = res[0], res[2]
                    dssp_predictions.append(dssp_code)

        end_time = time.time()
        dsspTime = end_time - start_time
        print(f"Time taken for DSSP prediction: {dsspTime:.2f} seconds.")

        with open(f"real_output.txt", "w") as dssp_out:
                print(f"DSSP Prediction Table \n", dssp_predictions, file=dssp_out)

        return dssp_predictions

    
    def _estimate_accuracy(self, custom_predictions, dssp_predictions):
        """
        Estimate the accuracy of secondary structure prediction compared to DSSP results.
        """

        matched = 0
        unknown = 0
        total_residues = len(dssp_predictions)
        
        for (i, custom_pred), dssp_pred in zip(custom_predictions, dssp_predictions):
            if dssp_pred == custom_pred or (dssp_pred == '-' and custom_pred == '-'):
                matched += 1
            elif dssp_pred == '-':
                unknown += 1
        
        accuracy_percentage = (matched / (matched+unknown)) * 100 if total_residues > 0 else 0
        print(f"Accuracy Percentage: {accuracy_percentage}%")
        print(f"Unknown Predictions Count: {unknown}")

        return accuracy_percentage


    def display_results_in_gui(self, predictions):
        """
        Displays secondary structure results using tkinter.
        """

        root = tk.Tk()
        root.title("Secondary Structure Prediction")

        columns = ('Residue', 'Structure')

        tree = ttk.Treeview(root, columns=columns, show='headings')
        tree.heading('Residue', text='Residue')
        tree.heading('Structure', text='Structure')

        for pred in predictions:
            tree.insert('', tk.END, values=pred)

        tree.pack(expand=True, fill=tk.BOTH)
        root.mainloop()


# Main
if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Use python script.py <PDB file path>")
        sys.exit(1)
    
    print("Thank you for trying out the DSSP Prototype for Protein Prediction (DPPP) !")
    
    # Calculate memory usage
    start_memory = memory_usage(-1, interval=0.1, timeout=1)

    # Specify PDB file
    pdb_file_path = sys.argv[1]

    dppp = Dppp(pdb_file_path)

    structure = dppp.read_pdb_file()

    dppp.get_structure_info(structure)

    # Get real DSSP predictions
    dssp_predictions = dppp.calculate_dssp_predictions(structure)

    # Get custom DSSP predictions
    custom_predictions = dppp.predict_multiprocessing(structure)

    # Compare and estimate accuracy
    dppp._estimate_accuracy(custom_predictions, dssp_predictions)

    # Show memory and time usage
    end_memory = memory_usage(-1, interval=0.1, timeout=1)
    print(f"Memory usage: {end_memory[0]- start_memory[0]} MiB")

    # Display results
    dppp.display_results_in_gui(custom_predictions)


