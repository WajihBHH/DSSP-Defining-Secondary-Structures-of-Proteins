

#!/usr/bin/env python3
"""
DSSP-Style Secondary Structure Assignment Script
Implements hydrogen bond-based assignment of alpha-helices and beta-sheets
"""
import math
import sys
import warnings
import numpy as np

from Bio import PDB
from Bio.PDB import PDBParser, DSSP
from collections import defaultdict



warnings.filterwarnings('ignore')

class SecondaryStructureAssigner:
    """
    DSSP-style secondary structure assignment based on hydrogen bonding patterns
    """
    
    def __init__(self, pdb_file):
        """
        Initialize the secondary structure assigner
        
        Args:
            pdb_file (str): Path to PDB file
        """
        self.pdb_file = pdb_file
        self.structure = None
        self.residues = []
        self.backbone_atoms = {}
        self.hydrogen_bonds = []
        self.assignments = {}
        
    def parse_pdb(self):
        """
        Parse PDB file and extract backbone atoms for each residue
        """
        print("Step 1: Parsing PDB file...")
        
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure('protein', self.pdb_file)
        
        # Extract residues and backbone atoms
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if residue.has_id('N') and residue.has_id('CA') and residue.has_id('C') and residue.has_id('O'):
                        res_id = (chain.id, residue.id[1])  # (chain, residue_number)
                        self.residues.append(res_id)
                        
                        # Store backbone atoms
                        self.backbone_atoms[res_id] = {
                            'N': residue['N'].get_coord(),
                            'CA': residue['CA'].get_coord(),
                            'C': residue['C'].get_coord(),
                            'O': residue['O'].get_coord()
                        }
        
        print(f"Parsed {len(self.residues)} residues with complete backbone atoms")
        
    def estimate_hydrogen_positions(self):
        """
        Estimate hydrogen atom positions for amide groups (N-H)
        N-H bond length : ~1.0 Å
        """
        print("Step 2: Estimating hydrogen atom positions...")
        
        for i, res_id in enumerate(self.residues):
            atoms = self.backbone_atoms[res_id]
            
            # For first residue, skip (no previous C to reference)
            if i == 0:
                # Place H along N-CA direction (approximate)
                n_pos = atoms['N']
                ca_pos = atoms['CA']
                direction = n_pos - ca_pos
                direction = direction / np.linalg.norm(direction)
                h_pos = n_pos + direction * 1.0  # 1.0 Å N-H bond
                
            else:
                # Use previous residue's C to determine H position
                prev_res_id = self.residues[i-1]
                prev_c = self.backbone_atoms[prev_res_id]['C']
                n_pos = atoms['N']
                ca_pos = atoms['CA']
                
                # H position: bisector of C(i-1)-N-CA angle, 1.0 Å from N
                v1 = prev_c - n_pos
                v2 = ca_pos - n_pos
                v1 = v1 / np.linalg.norm(v1)
                v2 = v2 / np.linalg.norm(v2)
                
                # Bisector direction (opposite to average of normalized vectors)
                bisector = -(v1 + v2)
                if np.linalg.norm(bisector) > 0:
                    bisector = bisector / np.linalg.norm(bisector)
                    h_pos = n_pos + bisector * 1.0
                else:
                    # Fallback: place along N-CA direction
                    direction = n_pos - ca_pos
                    direction = direction / np.linalg.norm(direction)
                    h_pos = n_pos + direction * 1.0
            
            # Store estimated H position
            self.backbone_atoms[res_id]['H'] = h_pos
            
        print("Estimated hydrogen positions for all residues")
    
    def calculate_hydrogen_bonds(self):
        """
        Calculate hydrogen bonds using Kabsch & Sander energy formula :
        E = 0.084 * (1/r_ON + 1/r_CH - 1/r_OH - 1/r_CN) * 332
        """
        print("Step 3: Calculating hydrogen bonds...")
        
        self.hydrogen_bonds = []
        energy_cutoff = -0.5  # kcal/mol
        
        for i, donor_id in enumerate(self.residues):
            for j, acceptor_id in enumerate(self.residues):
                # Skip same residue and adjacent residues
                if abs(i - j) < 2:
                    continue
                
                donor_atoms = self.backbone_atoms[donor_id]
                acceptor_atoms = self.backbone_atoms[acceptor_id]
                
                # Donor: N-H, Acceptor: C=O
                n_pos = donor_atoms['N']
                h_pos = donor_atoms['H']
                c_pos = acceptor_atoms['C']
                o_pos = acceptor_atoms['O']
                
                # Calculate distances (in Angstroms)
                r_on = np.linalg.norm(o_pos - n_pos)
                r_ch = np.linalg.norm(c_pos - h_pos)
                r_oh = np.linalg.norm(o_pos - h_pos)
                r_cn = np.linalg.norm(c_pos - n_pos)
                
                # Avoid division by zero
                if min(r_on, r_ch, r_oh, r_cn) < 0.1:
                    continue
                
                # Kabsch & Sander formula
                energy = 0.084 * (1/r_on + 1/r_ch - 1/r_oh - 1/r_cn) * 332

                if energy < energy_cutoff and r_oh < 3.5:  # Additional distance cutoff
                    self.hydrogen_bonds.append({
                        'donor': donor_id,
                        'acceptor': acceptor_id,
                        'donor_idx': i,
                        'acceptor_idx': j,
                        'energy': energy,
                        'distance': r_oh
                    })
        
        print(f"Found {len(self.hydrogen_bonds)} hydrogen bonds")
    
    def assign_secondary_structure(self):
        """
        Assign secondary structure based on hydrogen bonding patterns
        H = Alpha helix, E = Beta sheet, C = Coil
        """
        print("Step 4: Assigning secondary structures...")
        
        for res_id in self.residues:
            self.assignments[res_id] = 'C'
        
        hbond_map = defaultdict(list)
        for hb in self.hydrogen_bonds:
            hbond_map[hb['donor_idx']].append(hb['acceptor_idx'])
            hbond_map[hb['acceptor_idx']].append(hb['donor_idx'])
        
        # 1. Identify alpha-helices (i to i+4 hydrogen bonds)
        helix_residues = set()
        for hb in self.hydrogen_bonds:
            donor_idx = hb['donor_idx']
            acceptor_idx = hb['acceptor_idx']
            
            # Alpha-helix: N-H(i) ... O=C(i-4) or N-H(i) ... O=C(i+4)
            if abs(donor_idx - acceptor_idx) == 4:
                # Mark both residues and intermediate residues as helix
                start_idx = min(donor_idx, acceptor_idx)
                end_idx = max(donor_idx, acceptor_idx)
                for idx in range(start_idx, end_idx + 1):
                    helix_residues.add(idx)
        
        # Assign helix
        for idx in helix_residues:
            if idx < len(self.residues):
                res_id = self.residues[idx]
                self.assignments[res_id] = 'H'
        
        # 2. Identify beta-sheets (cross-strand hydrogen bonding)
        sheet_pairs = []
        for hb in self.hydrogen_bonds:
            donor_idx = hb['donor_idx']
            acceptor_idx = hb['acceptor_idx']
            
            # Beta-sheet: non-local hydrogen bonds (not helix pattern)
            if abs(donor_idx - acceptor_idx) > 4:
                sheet_pairs.append((donor_idx, acceptor_idx))
        
        # Find beta-sheet regions (simplified approach)
        sheet_residues = set()
        for donor_idx, acceptor_idx in sheet_pairs:
            found_partner = False
            for other_donor, other_acceptor in sheet_pairs:
                if (abs(donor_idx - other_donor) <= 2 and 
                    abs(acceptor_idx - other_acceptor) <= 2 and
                    (donor_idx != other_donor or acceptor_idx != other_acceptor)):
                    found_partner = True
                    break
            
            if found_partner:
                sheet_residues.add(donor_idx)
                sheet_residues.add(acceptor_idx)
        
        # Assign sheets (only if not already helix)
        for idx in sheet_residues:
            if idx < len(self.residues):
                res_id = self.residues[idx]
                if self.assignments[res_id] == 'C':  # Don't override helix
                    self.assignments[res_id] = 'E'
        
        # Count assignments
        counts = {'H': 0, 'E': 0, 'C': 0}
        for assignment in self.assignments.values():
            counts[assignment] += 1
        
        print(f"Assigned: {counts['H']} helix, {counts['E']} sheet, {counts['C']} coil")
    
    def get_reference_structure(self):
        """
        Get reference secondary structure from PDB HELIX/SHEET records or DSSP
        """
        print("Step 5: Getting reference secondary structure...")
        
        reference = {}
        
        # Initialize all as coil
        for res_id in self.residues:
            reference[res_id] = 'C'
        
        try:
            # Try to use DSSP if available
            dssp = DSSP(self.structure[0], self.pdb_file)
            
            for key in dssp.keys():
                chain_id = key[0]
                res_num = key[1][1]
                res_id = (chain_id, res_num)
                
                if res_id in self.residues:
                    dssp_code = dssp[key][2]
                    
                    # Convert DSSP codes to simplified H/E/C
                    if dssp_code in ['H', 'G', 'I']:  # Alpha-helix, 3-10 helix, pi-helix
                        reference[res_id] = 'H'
                    elif dssp_code in ['E', 'B']:  # Beta-strand, isolated beta-bridge
                        reference[res_id] = 'E'
                    else:
                        reference[res_id] = 'C'
            
            print("Used DSSP for reference structure")
            
        except:
            print("DSSP not available, parsing PDB records...")
            
            with open(self.pdb_file, 'r') as f:
                lines = f.readlines()
            
            # Parse HELIX records
            for line in lines:
                if line.startswith('HELIX'):
                    chain_id = line[19].strip()
                    start_res = int(line[21:25].strip())
                    end_res = int(line[33:37].strip())
                    
                    for res_num in range(start_res, end_res + 1):
                        res_id = (chain_id, res_num)
                        if res_id in self.residues:
                            reference[res_id] = 'H'
            
            # Parse SHEET records
            for line in lines:
                if line.startswith('SHEET'):
                    chain_id = line[21].strip()
                    start_res = int(line[22:26].strip())
                    end_res = int(line[33:37].strip())
                    
                    for res_num in range(start_res, end_res + 1):
                        res_id = (chain_id, res_num)
                        if res_id in self.residues:
                            reference[res_id] = 'E'
        
        return reference
    
    def calculate_accuracy(self):
        """
        Calculate accuracy metrics comparing predicted vs reference structures
        """
        print("Step 6: Calculating accuracy metrics...")
        
        reference = self.get_reference_structure()
        
        # Overall accuracy
        total = len(self.residues)
        correct = sum(1 for res_id in self.residues 
                     if self.assignments[res_id] == reference[res_id])
        accuracy = correct / total if total > 0 else 0
        
        # Per-class metrics
        metrics = {}
        for struct_type in ['H', 'E', 'C']:
            # True positives, false positives, false negatives
            tp = sum(1 for res_id in self.residues 
                    if self.assignments[res_id] == struct_type and reference[res_id] == struct_type)
            fp = sum(1 for res_id in self.residues 
                    if self.assignments[res_id] == struct_type and reference[res_id] != struct_type)
            fn = sum(1 for res_id in self.residues 
                    if self.assignments[res_id] != struct_type and reference[res_id] == struct_type)
            
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
            
            metrics[struct_type] = {
                'precision': precision,
                'recall': recall,
                'f1': f1,
                'tp': tp,
                'fp': fp,
                'fn': fn
            }
        
        return accuracy, metrics, reference
    
    def run_analysis(self):
        """
        Run the complete secondary structure analysis pipeline
        """
        print("Starting DSSP-Style Secondary Structure Assignment")
        print("=" * 60)
        
        # Execute pipeline
        self.parse_pdb()
        self.estimate_hydrogen_positions()
        self.calculate_hydrogen_bonds()
        self.assign_secondary_structure()
        
        # Calculate accuracy
        accuracy, metrics, reference = self.calculate_accuracy()
        
        # Display results
        self.display_results(accuracy, metrics, reference)
        
        return self.assignments, accuracy, metrics
    
    def display_results(self, accuracy, metrics, reference):
        """
        Display analysis results in a formatted way
        """
        print("\n" + "=" * 60)
        print("RESULTS SUMMARY")
        print("=" * 60)
        
        print(f"\n Overall Accuracy: {accuracy:.3f} ({accuracy*100:.1f}%)")
        
        print(f"\n Per-Class Metrics:")
        print(f"{'Class':<10} {'Precision':<10} {'Recall':<10} {'F1-Score':<10}")
        print("-" * 50)
        
        class_names = {'H': 'Helix', 'E': 'Sheet', 'C': 'Coil'}
        for struct_type in ['H', 'E', 'C']:
            m = metrics[struct_type]
            print(f"{class_names[struct_type]:<10} {m['precision']:<10.3f} {m['recall']:<10.3f} {m['f1']:<10.3f}")
        
        print(f"\n Structure Assignments (First 20 residues):")
        print(f"{'Residue':<10} {'Predicted':<10} {'Reference':<10} {'Match':<10}")
        print("-" * 50)
        
        for i, res_id in enumerate(self.residues[:20]):
            chain, res_num = res_id
            pred = self.assignments[res_id]
            ref = reference[res_id]
            match = "✅" if pred == ref else "❌"
            print(f"{chain}{res_num:<9} {pred:<10} {ref:<10} {match}")
        
        if len(self.residues) > 20:
            print(f"... and {len(self.residues) - 20} more residues")
        
        print(f"\n Hydrogen Bonds Found: {len(self.hydrogen_bonds)}")
        if self.hydrogen_bonds:
            print("   Sample bonds (first 5):")
            for i, hb in enumerate(self.hydrogen_bonds[:5]):
                donor = self.residues[hb['donor_idx']]
                acceptor = self.residues[hb['acceptor_idx']]
                print(f"   {i+1}. {donor[0]}{donor[1]} → {acceptor[0]}{acceptor[1]} "
                      f"(E={hb['energy']:.2f} kcal/mol, d={hb['distance']:.2f} Å)")


def main():
    """
    Main function to run the secondary structure analysis
    """
    # Default PDB file (can be changed)
    pdb_file ="1zaa.pdb"  # Change this to your PDB file path

    if len(sys.argv) > 1:
        pdb_file = sys.argv[1]
    
    print(f" Analyzing protein structure from: {pdb_file}")
    
    try:
        # Create analyzer and run analysis
        analyzer = SecondaryStructureAssigner(pdb_file)
        assignments, accuracy, metrics = analyzer.run_analysis()
        
        # Optional: Save results to file
        output_file = pdb_file.replace('.pdb', '_secondary_structure.txt')
        with open(output_file, 'w') as f:
            f.write("Residue\tChain\tNumber\tPredicted\n")
            for res_id in analyzer.residues:
                chain, res_num = res_id
                assignment = assignments[res_id]
                f.write(f"{chain}{res_num}\t{chain}\t{res_num}\t{assignment}\n")
        
        print(f"\n Results saved to: {output_file}")
        
    except FileNotFoundError:
        print(f"Error: PDB file '{pdb_file}' not found!")
        print("Please make sure the file exists and try again.")
        print("\nUsage: python script.py [pdb_file]")
        
    except Exception as e:
        print(f"Error during analysis: {str(e)}")




if __name__ == "__main__":
    pdb_file = "src/1zaa.pdb"  # Make sure this file exists
    
    print(f"Analyzing protein structure from: {pdb_file}")
    
    try:
        analyzer = SecondaryStructureAssigner(pdb_file)
        assignments, accuracy, metrics = analyzer.run_analysis()
        
        output_file = pdb_file.replace('.pdb', '_secondary_structure.txt')
        with open(output_file, 'w') as f:
            f.write("Residue\tChain\tNumber\tPredicted\n")
            for res_id in analyzer.residues:
                chain, res_num = res_id
                assignment = assignments[res_id]
                f.write(f"{chain}{res_num}\t{chain}\t{res_num}\t{assignment}\n")
        
        print(f"\n Results saved to: {output_file}")
        
    except FileNotFoundError:
        print(f"Error: PDB file '{pdb_file}' not found!")
        
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
