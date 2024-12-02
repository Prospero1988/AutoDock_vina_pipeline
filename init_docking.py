#!/usr/bin/env python3
"""
This script automates the process of docking multiple ligands to multiple receptor proteins using AutoDock Vina.
It reads receptors PDB ID from a CSV file, downloads each receptor structure from the PDB database, prepares them for docking, and processes multiple ligands from an SDF file.

Technologies used:
- Biopython for handling PDB files and interacting with the PDB database.
- RDKit for cheminformatics operations.
- Open Babel for file format conversions and molecule preparation.
- PyMOL for molecular visualization and coordinate manipulation.
- P2Rank for binding site prediction.
- AutoDock Vina for molecular docking.

How to run:

python dock.py --pdb_ids receptors.csv --ligands ligand_file.sdf

Arguments:

- --pdb_ids: The name of the CSV file containing PDB IDs of receptor proteins, located in the ./receptors directory.
- --ligands: The name of the SDF file containing ligands to dock, located in the ./ligands directory.
- --tol: Tolerance in Ångströms to expand the docking pocket dimensions beyond those defined by P2Rank (default: 0).
- --pckt: Pocket number to use from P2Rank predictions (default: 1).
- --exhaust: Specifies how thorough the search should be for the best binding poses.
              Higher values increase precision but require more computation time (default: 20).
- --energy_range: Determines the range of energy scores (in kcal/mol) for poses to be considered (default: 2).

All receptor-related files will be saved in the ./PDB_ID directory. Each ligand's docking results will
be saved in ./PDB_ID/ligand_name_or_number.

Ensure that all required tools and libraries are installed and properly configured in your environment.
"""

import argparse
import os
import sys
import shutil
import subprocess
import logging
import re

from Bio import PDB
from Bio.PDB import PDBList
import pandas as pd
import numpy as np
import pymol
from pymol import cmd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Paths to external tools
P2RANK_PATH = "/usr/local/bin/prank"
VINA_PATH = "/usr/local/bin/vina_1.2.5_linux_x86_64"
OBABEL_PATH = "/usr/bin/obabel"

# Check if external tools are available
if not os.path.exists(P2RANK_PATH):
    raise FileNotFoundError(f"P2Rank executable not found at {P2RANK_PATH}")
if not os.path.exists(VINA_PATH):
    raise FileNotFoundError(f"AutoDock Vina executable not found at {VINA_PATH}")
if shutil.which(OBABEL_PATH) is None:
    raise EnvironmentError(f"Open Babel not found at {OBABEL_PATH}")

# Initialize PyMOL in quiet mode without GUI
pymol.finish_launching(['pymol', '-qc'])

# Logger setup
def logger_decorator(func):
    def wrapper(*args, **kwargs):
        logging.info(f"Running {func.__name__}")
        return func(*args, **kwargs)
    return wrapper

# Main function
def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Automated docking script for multiple ligands and receptors.")
    parser.add_argument('--pdb_ids', required=True, help='Name of the CSV file containing PDB IDs of receptor proteins, located in ./receptors.')
    parser.add_argument('--ligands', required=True, help='Name of the SDF file containing ligands, located in ./ligands.')
    parser.add_argument('--tol', type=int, default=0, help='Tolerance in Ångströms to expand the docking pocket dimensions beyond those defined by P2Rank (default: 0).')
    parser.add_argument('--pckt', type=int, default=1, help='Pocket number to use from P2Rank predictions (default: 1).')
    parser.add_argument('--exhaust', type=int, default=20, help='Specifies how thorough the search should be for the best binding poses. Higher values increase precision but require more computation time (default: 20).')
    parser.add_argument('--energy_range', type=int, default=2, help='Determines the range of energy scores (in kcal/mol) for poses to be considered (default: 2).')

    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_ids_file = args.pdb_ids
    receptors_file_path = os.path.join(script_dir, 'receptors', pdb_ids_file)

    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')

    # Folder for ligands
    ligands_folder = os.path.join(script_dir, 'ligands')
    ligands_file = args.ligands
    ligands_path = os.path.join(ligands_folder, ligands_file)

    if not os.path.exists(ligands_path):
        raise FileNotFoundError(f"Ligand file not found: {ligands_path}")
    else:
        logging.info(f"Ligand file found: {ligands_path}")

    suppl = Chem.SDMolSupplier(ligands_path)
    if not suppl:
        raise ValueError(f"Could not read ligands from {ligands_path} or file is empty.")
    else:
        logging.info(f"Number of ligands read: {len(suppl)}")

    # Read PDB IDs from the CSV file
    if not os.path.exists(receptors_file_path):
        raise FileNotFoundError(f"Receptors file not found: {receptors_file_path}")
    else:
        logging.info(f"Receptors file found: {receptors_file_path}")

    with open(receptors_file_path, 'r') as receptors_file:
        pdb_ids = [line.strip() for line in receptors_file if line.strip()]

    if not pdb_ids:
        raise ValueError("No PDB IDs found in the receptors file.")

    # Loop to work on multiple receptors
    for PDB_ID in pdb_ids:
        pocket_number = args.pckt
        receptor_name = PDB_ID
        tol = args.tol
        exhaustiveness = args.exhaust
        energy_range = args.energy_range

        # Set up paths
        receptor_folder = os.path.join(script_dir, receptor_name)
        if not os.path.exists(receptor_folder):
            os.makedirs(receptor_folder, exist_ok=True)

        # Set up logging for this receptor
        log_file = os.path.join(receptor_folder, f"{receptor_name}_console_output.log")
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s %(message)s')
        file_handler.setFormatter(formatter)
        logging.getLogger().addHandler(file_handler)

        logging.info(f"Processing receptor {PDB_ID}")
        logging.info(f"Tolerance (Å): {tol}")
        print(f"Docking tolerance set to {tol} Å.")

        # Output results file
        results_file = os.path.join(receptor_folder, f"{receptor_name}_results.txt")

        # Initialize list to collect ligand results
        ligand_results = []

        # Start processing
        try:
            # Download and prepare receptor
            downloaded_pdb_path, protein_name = download_pdb(PDB_ID, receptor_folder)
            dirty_pdb = os.path.join(receptor_folder, f'{receptor_name}_dirty.pdb')
            shutil.move(downloaded_pdb_path, dirty_pdb)

            fixed_pdb = os.path.join(receptor_folder, f'{receptor_name}_fixed.pdb')
            fix_pdb(dirty_pdb, fixed_pdb, ph=7.4)

            receptor_pdbqt = os.path.join(receptor_folder, f"{receptor_name}.pdbqt")
            prepare_receptor(fixed_pdb, receptor_pdbqt)

            # Run P2Rank
            output_dir = os.path.join(receptor_folder, '01_p2rank_output')
            run_p2rank(fixed_pdb, output_dir)

            # Get docking box parameters
            center_x, center_y, center_z, Size_x, Size_y, Size_z, predictions_csv = get_docking_box(output_dir, fixed_pdb, tol, pocket_number)

            with open(results_file, 'w') as rf:
                for idx, mol in enumerate(suppl):
                    if mol is None:
                        logging.warning(f"Skipping invalid molecule at index {idx}")
                        continue
                    ligand_name = mol.GetProp('_Name') if mol.HasProp('_Name') else f"ligand_{idx}"
                    ligand_folder = os.path.join(receptor_folder, ligand_name)
                    os.makedirs(ligand_folder, exist_ok=True)

                    ligand_pdb = os.path.join(ligand_folder, f"{ligand_name}.pdb")
                    ligand_pdbqt = os.path.join(ligand_folder, f"{ligand_name}.pdbqt")
                    output_pdbqt = os.path.join(ligand_folder, f"{receptor_name}_{ligand_name}.pdbqt")

                    # Save ligand to PDB
                    write_mol_to_pdb(mol, ligand_pdb)

                    # Generate 2D coordinates and save image
                    AllChem.Compute2DCoords(mol)
                    image_filename = os.path.join(ligand_folder, f"{ligand_name}.svg")
                    draw_molecule_to_file(mol, image_filename)

                    # Prepare ligand
                    prepare_ligand(ligand_pdb, ligand_pdbqt)

                    # Run docking
                    vina_output, affinities = run_vina(receptor_pdbqt, ligand_pdbqt, output_pdbqt,
                                        center_x, center_y, center_z,
                                        Size_x, Size_y, Size_z,
                                        exhaustiveness, energy_range)

                    # Save vina output to results file
                    rf.write(f"Ligand: {ligand_name}\n")
                    rf.write(vina_output)
                    rf.write("\n\n")

                    # Display information in the terminal after docking is complete
                    print(f'Ligand {ligand_name} docked successfully!')

                    # Generate visualization 
                    generate_visualizations(receptor_pdbqt, output_pdbqt, ligand_folder, receptor_name, ligand_name)
                    docking_image_path = os.path.join(ligand_folder, f"{receptor_name}_{ligand_name}_docking.png")

                    # Get the first affinity value
                    if affinities:
                        affinity = affinities[0][0]
                    else:
                        affinity = None

                    # Collect ligand results
                    ligand_results.append({
                        'name': ligand_name,
                        'image': image_filename,
                        'affinity': affinity,
                        'output_pdbqt': output_pdbqt,
                        'docking_image': docking_image_path
                    })

            # Generate HTML results file
            html_file = os.path.join(receptor_folder, f"{receptor_name}_results.html")
            generate_html_results(html_file, receptor_name, ligands_file, ligand_results, predictions_csv, protein_name, args.pckt, receptor_pdbqt)

            print(f'Results HTML file generated: {html_file}')

        except Exception as e:
            logging.error(f"An error occurred with PDB ID {PDB_ID}: {e}")
            print(f"An error occurred with PDB ID {PDB_ID}: {e}")
            # Remove the file handler after processing each receptor
            logging.getLogger().removeHandler(file_handler)
            continue  # Proceed to the next PDB_ID

        # Remove the file handler after processing each receptor
        logging.getLogger().removeHandler(file_handler)

# Function definitions (outside of the main function and loops)

@logger_decorator
def download_pdb(pdb_id, download_dir):
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)
    pdbl = PDBList()
    try:
        pdb_file_path = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb', pdir=download_dir)
    except Exception as e:
        logging.error(f"Failed to download PDB {pdb_id}: {e}")
        raise
    # Handle compressed files
    if pdb_file_path.endswith('.gz'):
        import gzip
        with gzip.open(pdb_file_path, 'rb') as f_in:
            uncompressed_path = pdb_file_path[:-3]  # Remove '.gz'
            with open(uncompressed_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(pdb_file_path)
        pdb_file_path = uncompressed_path
    # Rename to standard PDB file name
    expected_filename = os.path.join(download_dir, f"{pdb_id}.pdb")
    os.rename(pdb_file_path, expected_filename)
    if not os.path.exists(expected_filename):
        raise FileNotFoundError(f"Failed to download PDB {pdb_id}.")
    logging.info(f"PDB file downloaded: {expected_filename}")
    
    # Extract protein name
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, expected_filename)
    header = structure.header  # Header contains metadata
    protein_name = header.get('name', 'Unknown')  # Extract protein name if available
    logging.info(f"Protein name for {pdb_id}: {protein_name}")
    
    return expected_filename, protein_name

@logger_decorator
def fix_pdb(input_pdb, output_pdb, ph=7.4):
    try:
        fixer = PDBFixer(filename=input_pdb)

        # Get all chains and their residue counts using chain indices
        chain_residue_counts = {}
        chain_list = list(fixer.topology.chains())
        for index, chain in enumerate(chain_list):
            residue_count = sum(1 for residue in chain.residues())
            chain_residue_counts[index] = residue_count
            chain_id = chain.id
            print(f"Found chain index: {index}, id: '{chain_id}', residues: {residue_count}")

        if not chain_residue_counts:
            raise ValueError("No chains found in the PDB structure.")

        # Decide which chain to keep (chain with the highest number of residues)
        chain_to_keep_index = max(chain_residue_counts, key=chain_residue_counts.get)
        chain_to_keep_id = chain_list[chain_to_keep_index].id
        logging.info(f"Keeping chain index {chain_to_keep_index} (id '{chain_to_keep_id}') with {chain_residue_counts[chain_to_keep_index]} residues.")
        print(f"Keeping chain index {chain_to_keep_index} (id '{chain_to_keep_id}') with {chain_residue_counts[chain_to_keep_index]} residues.")

        # Remove all other chains by their indices
        chains_to_remove = [i for i in range(len(chain_list)) if i != chain_to_keep_index]
        fixer.removeChains(chains_to_remove)

        # Remove heterogens (including water)
        fixer.removeHeterogens(keepWater=False)

        # Find missing residues
        fixer.findMissingResidues()

        # Find missing atoms and add them
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        # Add missing hydrogens
        fixer.addMissingHydrogens(ph)

        # Save the fixed PDB file
        with open(output_pdb, 'w') as outfile:
            PDBFile.writeFile(fixer.topology, fixer.positions, outfile)
        logging.info(f"Fixed PDB saved as {output_pdb}")
    except Exception as e:
        logging.error(f"Error in fixing PDB: {e}")
        print(f"Error in fixing PDB: {e}")
        raise

@logger_decorator
def prepare_receptor(input_pdb, output_pdbqt):
    try:
        subprocess.run([OBABEL_PATH, "-i", "pdb", input_pdb, "-o", "pdbqt", "-O", output_pdbqt, "-xr"], check=True)
        logging.info(f"Receptor prepared: {output_pdbqt}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in preparing receptor: {e}")
        print(f"Error in preparing receptor: {e}")
        raise

@logger_decorator
def prepare_ligand(input_pdb, output_pdbqt):
    try:
        subprocess.run([OBABEL_PATH, "-i", "pdb", input_pdb, "-o", "pdbqt", "-O", output_pdbqt, "-h"], check=True)
        logging.info(f"Ligand prepared: {output_pdbqt}")
        print(f"Ligand {os.path.basename(input_pdb)} prepared successfully!")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in preparing ligand: {e}")
        print(f"Error in preparing ligand: {e}")
        raise

@logger_decorator
def run_p2rank(receptor_pdb, output_dir):
    try:
        subprocess.run([P2RANK_PATH, 'predict', '-f', receptor_pdb, '-o', output_dir], check=True)
        logging.info("P2Rank prediction completed.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in P2Rank prediction: {e}")
        print(f"Error in P2Rank prediction: {e}")
        raise

@logger_decorator
def get_docking_box(output_dir, receptor_pdb, tol, pocket_number):
    try:
        receptor_base_name = os.path.basename(receptor_pdb).split('.')[0]
        predictions_csv = os.path.join(output_dir, f'{receptor_base_name}.pdb_predictions.csv')
        residues_csv = os.path.join(output_dir, f'{receptor_base_name}.pdb_residues.csv')

        df = pd.read_csv(predictions_csv)
        pred = pd.read_csv(residues_csv)

        # Remove whitespace characters from column names
        df.columns = df.columns.str.strip()
        pred.columns = pred.columns.str.strip()

        # Obtaining center_x, center_y, center_z values
        center_x = float(df['center_x'].iloc[pocket_number - 1])
        center_y = float(df['center_y'].iloc[pocket_number - 1])
        center_z = float(df['center_z'].iloc[pocket_number - 1])

        # Selection of residues for chosen pocket 
        pocket = pred[pred['pocket'] == pocket_number]

        if pocket.empty:
            raise ValueError(f"P2Rank did not predict pocket number {pocket_number}")

        logging.info(f"Using pocket number {pocket_number} from P2Rank predictions.")
        print(f"Using pocket number {pocket_number} from P2Rank predictions.")

        resi_numbers = [str(res).strip() for res in pocket['residue_label']]
        resi = '+'.join(resi_numbers)

        cmd.load(receptor_pdb)
        cmd.select('pocket', f'resi {resi}')
        cmd.show('cartoon')

        alpha_carbons = []
        model = cmd.get_model('pocket and name CA')
        for atom in model.atom:
            alpha_carbons.append([atom.coord[0], atom.coord[1], atom.coord[2]])

        alpha_carbons = np.array(alpha_carbons)
        min_coords = np.min(alpha_carbons, axis=0)
        max_coords = np.max(alpha_carbons, axis=0)
        cube_size = max_coords - min_coords
        Size_x, Size_y, Size_z = cube_size + np.array([tol, tol, tol])

        logging.info(f"Docking box dimensions before tolerance: {cube_size}")
        logging.info(f"Docking box dimensions after tolerance: {Size_x}, {Size_y}, {Size_z}")
        logging.info(f"Tolerance applied: {tol}")
        logging.info(f"Docking box parameters obtained.")
        cmd.delete('all')
        return center_x, center_y, center_z, Size_x, Size_y, Size_z, predictions_csv
    except Exception as e:
        logging.error(f"Error in getting docking box parameters: {e}")
        print(f"Error in getting docking box parameters: {e}")
        raise

@logger_decorator
def run_vina(receptor_pdbqt, ligand_pdbqt, output_pdbqt,
             center_x, center_y, center_z,
             size_x, size_y, size_z,
             exhaustiveness, energy_range):
    try:
        vina_command = [
            VINA_PATH,
            '--receptor', receptor_pdbqt,
            '--ligand', ligand_pdbqt,
            '--out', output_pdbqt,
            '--energy_range', str(energy_range),
            '--center_x', str(center_x),
            '--center_y', str(center_y),
            '--center_z', str(center_z),
            '--size_x', str(size_x),
            '--size_y', str(size_y),
            '--size_z', str(size_z),
            '--exhaustiveness', str(exhaustiveness),
            '--seed', '1988',  # Optionally, to have reproducible results

        ]

        # Call Vina without capturing the output
        subprocess.run(vina_command, check=True)
        logging.info(f"Docking completed for ligand: {ligand_pdbqt}")

        # Read the docking results from the output_pdbqt file
        affinities = []
        with open(output_pdbqt, 'r') as f:
            for line in f:
                if line.startswith('REMARK VINA RESULT:'):
                    parts = line.strip().split()
                    if len(parts) >= 6:
                        # parts[3]: affinity
                        # parts[4]: RMSD l.b.
                        # parts[5]: RMSD u.b.
                        affinity = float(parts[3])
                        rmsd_lb = float(parts[4])
                        rmsd_ub = float(parts[5])
                        affinities.append((affinity, rmsd_lb, rmsd_ub))
        if affinities:
            # Create a table similar to Vina's output
            table_output = 'mode |   affinity | rmsd l.b.| rmsd u.b.\n'
            table_output += '-----+------------+----------+----------\n'
            for idx, (affinity, rmsd_lb, rmsd_ub) in enumerate(affinities, start=1):
                table_output += f'{idx:>4}    {affinity:>10.4f}   {rmsd_lb:>8.3f}   {rmsd_ub:>8.3f}\n'
        else:
            table_output = "No docking results found."

        return table_output, affinities  # Modified to return list of affinities

    except subprocess.CalledProcessError as e:
        logging.error(f"Error in docking with Vina: {e}")
        print(f"Error in docking with Vina: {e}")
        raise

@logger_decorator
def generate_visualizations(receptor_pdbqt, output_pdbqt, output_folder, receptor_name, ligand_name):
    try:
        # Load structures into PyMOL from PDBQT files (containing charges)
        cmd.load(receptor_pdbqt, 'receptor')
        cmd.load(output_pdbqt, 'ligand')
        
        # Add hydrogen atoms (if not present)
        cmd.h_add('receptor')
        cmd.h_add('ligand')
        
        # Color the receptor and ligand
        cmd.color('cyan', 'receptor')
        cmd.color('red', 'ligand')
    
        # Set representations
        cmd.show('cartoon', 'receptor')
        cmd.show('sticks', 'ligand')
    
        # Hide hydrogen atoms in the visualization
        cmd.hide('sticks', '(elem H)')
    
        # Define donors and acceptors in the receptor and ligand
        cmd.select('donors_receptor', '(receptor) and (elem N,O) and (neighbor hydro)')
        cmd.select('acceptors_receptor', '(receptor) and (elem N,O) and (not neighbor hydro)')
        cmd.select('donors_ligand', '(ligand) and (elem N,O) and (neighbor hydro)')
        cmd.select('acceptors_ligand', '(ligand) and (elem N,O) and (not neighbor hydro)')
    
        # Set cutoff distance for hydrogen bonds (in Ångströms)
        hbond_cutoff = 3.5  # You can adjust this value as needed
    
        # Find and display hydrogen bonds between the ligand and receptor
        # Donors from ligand to acceptors in receptor
        cmd.distance('hbonds', 'donors_ligand', 'acceptors_receptor', hbond_cutoff)
    
        # Donors from receptor to acceptors in ligand
        cmd.distance('hbonds', 'donors_receptor', 'acceptors_ligand', hbond_cutoff)
    
        # Hide distance labels
        cmd.set('label_size', 0)
        cmd.hide('labels')
    
        # **Add visualization of salt bridges using charges from PDBQT**
    
        # Define positively charged residues in the receptor (e.g., Lys, Arg)
        cmd.select('positive_residues_receptor', 'receptor and ((resn LYS and name NZ) or (resn ARG and name NH1+NH2+NE))')
    
        # Define negatively charged residues in the receptor (e.g., Asp, Glu)
        cmd.select('negative_residues_receptor', 'receptor and ((resn ASP and name OD1+OD2) or (resn GLU and name OE1+OE2))')
    
        # Define charged atoms in the ligand based on charges from PDBQT
        cmd.select('positive_atoms_ligand', 'ligand and partial_charge > 0.2')
        cmd.select('negative_atoms_ligand', 'ligand and partial_charge < -0.2')
    
        # Set cutoff distance for salt bridges (in Ångströms)
        salt_bridge_cutoff = 4.0  # You can adjust this value as needed
    
        # Find and display salt bridges between the ligand and receptor
        # Positively charged atoms from ligand to negatively charged residues in receptor
        cmd.distance('salt_bridges', 'positive_atoms_ligand', 'negative_residues_receptor', salt_bridge_cutoff)
    
        # Negatively charged atoms from ligand to positively charged residues in receptor
        cmd.distance('salt_bridges', 'negative_atoms_ligand', 'positive_residues_receptor', salt_bridge_cutoff)
    
        # Set color and style for salt bridges
        cmd.set_color('salt_bridge_color', [1.0, 0.5, 0.0])  # Orange color
        cmd.color('salt_bridge_color', 'salt_bridges')
        cmd.set('dash_width', 4.0, 'salt_bridges')  # Thicker lines
        cmd.set('dash_gap', 0.0, 'salt_bridges')    # Continuous lines
    
        # Adjust visualization and set the camera
        set_visualization_and_focus()
    
        # Save the image as PNG
        image_path = os.path.join(output_folder, f"{receptor_name}_{ligand_name}_docking.png")
        cmd.ray(1920, 1080)
        cmd.png(image_path, width=1920, height=1080, dpi=300)
        logging.info(f"Visualization saved: {image_path}")
    
        # Clear the PyMOL session
        cmd.delete('all')
    except Exception as e:
        logging.error(f"Error in generating visualization: {e}")
        print(f"Error in generating visualization: {e}")
        raise


def set_visualization_and_focus():
    cmd.show('cartoon', 'all')
    cmd.bg_color('white')
    cmd.center('ligand')  # Centring on the ligand
    cmd.zoom('ligand', buffer=1.0)  # Zoom in on the ligand 
    cmd.move('z', -50)  # Odsuń kamerę o 10 jednostek w osi Z

    # Deleting clipping settings
    # cmd.clip('near', -50)
    # cmd.clip('far', 50)
    cmd.set('ray_trace_mode', 1)
    cmd.set('orthoscopic', 1)
    cmd.set('ray_trace_frames', 1)
    cmd.set('antialias', 2)
    cmd.set('ambient', 0.5)
    cmd.set('spec_reflect', 0.5)
    cmd.set('specular', 1)


def write_mol_to_pdb(mol, pdb_filename):
    try:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        Chem.MolToPDBFile(mol, pdb_filename)
        logging.info(f"Ligand PDB saved: {pdb_filename}")
    except Exception as e:
        logging.error(f"Error in writing molecule to PDB: {e}")
        print(f"Error in writing molecule to PDB: {e}")
        raise

def draw_molecule_to_file(mol, image_filename):
    try:
        # Prepare SVG drawer
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 150)
        options = drawer.drawOptions()
        options.padding = 0.1  # Small margin
        options.fixedFontSize = 11  # Fixed font size
        options.useFixedFontSize = True  # Enforce fixed font size
        options.minFontSize = 6  # Minimum font size
        options.bondLineWidth = 2  # Thicker bond lines
        drawer.SetDrawOptions(options)

        # Draw molecule
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()

        # Remove XML declaration and fix namespace
        svg = svg.replace('<?xml version=\'1.0\' encoding=\'utf-8\'?>\n', '')
        svg = svg.replace('xmlns:svg=', 'xmlns=')

        # Save SVG file
        with open(image_filename, 'w') as f:
            f.write(svg)
        logging.info(f"Ligand image saved: {image_filename}")
    except Exception as e:
        logging.error(f"Error in drawing molecule image: {e}")
        print(f"Error in drawing molecule image: {e}")
        raise

def generate_html_results(html_file, receptor_name, ligands_file, ligand_results, predictions_csv, protein_name, pckt, receptor_pdbqt):
    try:
        # Load data from P2RANK and clean whitespace
        p2rank_csv = predictions_csv
        df_p2rank = pd.read_csv(p2rank_csv)
        df_p2rank.columns = df_p2rank.columns.str.strip()
        df_p2rank = df_p2rank.applymap(lambda x: x.strip() if isinstance(x, str) else x)
        df_p2rank['score'] = df_p2rank['score'].astype(float).map("{:.2f}".format)
        df_p2rank['probability'] = df_p2rank['probability'].astype(float).map("{:.2f}".format)
    
        with open(html_file, 'w', encoding='utf-8') as hf:
            hf.write('<html>\n')
            hf.write('<head>\n')
            hf.write('<title>Docking Results</title>\n')
            hf.write('<style>\n')
            hf.write('body { background-color: white; font-family: Arial, sans-serif; }\n')
            hf.write('table { border-collapse: collapse; margin: auto; }\n')
            hf.write('th, td { border: 1px solid black; padding: 5px; text-align: center; vertical-align: middle; }\n')
            hf.write('th { background-color: #f2f2f2; }\n')
            hf.write('img { display: block; margin: auto; }\n')
            hf.write('.probability { background-color: #ffecd9; } /* Pastel orange */\n')
            hf.write('.docking-energy { background-color: #dfffe0; } /* Pastel green */\n')  # Add pastel green styling
            # Style for the "Docking Results" column in the first table
            hf.write('td:nth-child(6), th:nth-child(6) {\n')
            hf.write('  max-width: 200px;\n')  # Maximum width
            hf.write('  word-wrap: break-word;\n')  # Word wrapping
            hf.write('  white-space: normal;\n')  # Normal white space
            hf.write('  text-align: left;\n')  # Align left
            hf.write('}\n')
            # Style for the "score" column in the second table with additional margins
            hf.write('.p2rank-table td:nth-child(3), .p2rank-table th:nth-child(3) {\n')
            hf.write('  white-space: nowrap;\n')  # Adjust width to content
            hf.write('  padding-left: 15px;\n')  # Left padding
            hf.write('  padding-right: 15px;\n')  # Right padding
            hf.write('}\n')
            # Style for the "residue_ids" column in the second table
            hf.write('.p2rank-table td:nth-child(5), .p2rank-table th:nth-child(5) {\n')
            hf.write('  max-width: 400px;\n')  # Ograniczenie szerokości kolumny
            hf.write('  word-wrap: break-word;\n')  # Zawijanie tekstu w kolumnie
            hf.write('  white-space: normal;\n')  # Normalne białe znaki dla zawijania
            hf.write('  padding: 15px;\n')  # Wewnętrzne marginesy
            hf.write('  text-align: left;\n')  # Opcjonalnie: wyrównanie tekstu do lewej
            hf.write('}\n')
            hf.write('</style>\n')
    
            hf.write('</head>\n')
            hf.write('<body>\n')
    
            # Header for the docking results table
            receptor_pdbqt = os.path.basename(receptor_pdbqt)
            header_text = f'Docking results for receptor with PDB code: <span style="color: red;">{receptor_name}</span></br>using structures from file: <span style="color: navy;">{ligands_file}</span></br>'
            hf.write(f'<h2 style="text-align: center;">{header_text}</h2>\n')
            hf.write(f'</br><h2 style="text-align: center; color: green; max-width: 600px; word-wrap: break-word; margin: auto;">{protein_name}</h2></br>')
            hf.write(f'<h3 style="text-align: center;">Docking to pocket number: <span style="color: red;">{pckt}</span></br></h2>\n')
            hf.write(f'<div style="text-align: center;">')
            hf.write(f'<a href="{receptor_pdbqt}" download="receptor_structure.pdbqt" type="application/octet-stream">Receptor structure (file .PDBQT). DOWNLOAD</a>')
            hf.write('</div>')
            hf.write('</br>')
            
            # First table: Docking results
            hf.write('<table>\n')
            hf.write('<tr><th>Number</th><th>Compound Name</th><th>Structure</th><th>Docking Image</th><th class="docking-energy">Docking Energy<br/>(kcal/mol)</th><th>Docking Results</th></tr>\n')
            for idx, result in enumerate(ligand_results, start=1):
                name = result['name']
                image_path = os.path.relpath(result['image'], os.path.dirname(html_file))
                docking_image_path = os.path.relpath(result['docking_image'], os.path.dirname(html_file))
                affinity = result['affinity']
                if affinity is not None:
                    affinity_str = f"{affinity:.2f}"
                else:
                    affinity_str = 'N/A'
                # Path to the output file
                output_pdbqt_path = os.path.relpath(result['output_pdbqt'], os.path.dirname(html_file))
                link_text = os.path.basename(result['output_pdbqt'])
                hf.write('<tr>\n')
                hf.write(f'<td>{idx}</td>\n')  # 1st column
                hf.write(f'<td>{name}</td>\n')  # 2nd column
                hf.write(f'<td><img src="{image_path}" alt="{name}" width="400"/></td>\n')  # 3rd column
                # New column for docking image (4th column)
                hf.write(f'<td><a href="{docking_image_path}" target="_blank"><img src="{docking_image_path}" alt="Docking Image" style="max-width: 150px; max-height: 150px;"></a></td>\n')
                hf.write(f'<td class="docking-energy">{affinity_str}</td>\n')  # 5th column
                hf.write(f'<td><a href="{output_pdbqt_path}" download="{link_text}" type="application/octet-stream">{link_text}</a></td>\n')  # 6th column
                hf.write('</tr>\n')
            hf.write('</table>\n')
            hf.write('</br>')
    
            # Link to detailed results
            results_file = f"{receptor_name}_results.txt"
            hf.write('<p style="text-align: center; margin-top: 20px;">\n')
            hf.write(f'<a href="{results_file}" target="_blank">Detailed results for each compound (All docking poses). CLICK</a>\n')
            hf.write('</p>\n')
            hf.write('</br></br>')
    
            # Header for the table with P2RANK data
            p2rank_header = f'P2RANK: identified docking pockets for receptor with PDB code: {receptor_name}'
            hf.write(f'<h3 style="text-align: center; margin-top: 20px;">{p2rank_header}</h3>\n')
    
            # Second table: P2RANK data
            hf.write('<table class="p2rank-table">\n')  # Added class "p2rank-table"
            hf.write('<tr>')
            for col in ['name', 'rank', 'score', 'probability', 'residue_ids']:
                hf.write(f'<th>{col}</th>')
            hf.write('</tr>\n')
            for _, row in df_p2rank.iterrows():
                hf.write('<tr>')
                for col in ['name', 'rank', 'score', 'probability', 'residue_ids']:
                    value = row[col]
                    if col == 'probability':
                        hf.write(f'<td class="probability">{value}</td>')
                    else:
                        hf.write(f'<td>{value}</td>')
                hf.write('</tr>\n')
            hf.write('</table>\n')
            hf.write('</br>')
            # Link to detailed results
            residues_csv_file = os.path.join('01_p2rank_output', f'{receptor_name}_fixed.pdb_residues.csv')
            hf.write('<p style="text-align: center; margin-top: 20px;">\n')
            hf.write(f'<a href="{residues_csv_file}" target="_blank">Detailed information about individual amino acids </br>and their involvement in docking pockets. CLICK</a>\n')
            hf.write('</p>\n')
    
            hf.write('</br></br>')
    
            # Author information at the end
            hf.write('<p style="font-size: small; text-align: center; margin-top: 20px;">\n')
            hf.write('Docking system based on <b>AutoDock Vina v.1.2.5</b> and <b>P2RANK v.2.4.2</b><br/>\n')
            hf.write('<b>Author:</b> Arkadiusz Leniak <b>email:</b> arkadiusz.leniak@gmail.com<br/>\n')
            hf.write('<b>github:</b> <a href="https://github.com/Prospero1988/AutoDock_vina_pipeline">https://github.com/Prospero1988/AutoDock_vina_pipeline</a>\n')
            hf.write('</p>\n')
    
            hf.write('</body>\n')
            hf.write('</html>\n')
        logging.info(f"HTML results file generated: {html_file}")
    except Exception as e:
        logging.error(f"Error in generating HTML results: {e}")
        print(f"Error in generating HTML results: {e}")
        raise


if __name__ == "__main__":
    main()
