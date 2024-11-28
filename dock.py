#!/usr/bin/env python3
"""
This script automates the process of docking multiple ligands to a receptor protein using AutoDock Vina.
It downloads the receptor structure from the PDB database, prepares it for docking, and processes multiple
ligands from an SDF file. The script utilizes various bioinformatics tools and libraries such as Biopython,
RDKit, Open Babel, PyMOL, and P2Rank.

Technologies used:
- Biopython for handling PDB files and interacting with the PDB database.
- RDKit for cheminformatics operations.
- Open Babel for file format conversions and molecule preparation.
- PyMOL for molecular visualization and coordinate manipulation.
- P2Rank for binding site prediction.
- AutoDock Vina for molecular docking.

How to run:

python dock.py --pdb_id PDB_ID --ligands ligand_file.sdf

Arguments:

- --pdb_id: The PDB ID of the receptor protein to download and use for docking.
- --ligands: The name of the SDF file containing ligands to dock, located in the ./ligands directory.
- --tol: Toleration in Ångströms to expand the docking pocket dimensions beyond those defined by P2RANK. (default: 0)
- --pckt: Pocket number to use from P2Rank predictions. (default: 1)


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
    parser = argparse.ArgumentParser(description="Automated docking script for multiple ligands.")
    parser.add_argument('--pdb_id', required=True, help='PDB ID of the receptor protein.')
    parser.add_argument('--ligands', required=True, help='Name of the SDF file containing ligands, located in ./ligands.')
    parser.add_argument('--tol', type=int, default=0, help='Toleration in Ångströms to expand the docking pocket dimensions beyond those defined by P2RANK.')
    parser.add_argument('--pckt', type=int, default=1, help='Pocket number to use from P2Rank predictions (default: 1)')

    args = parser.parse_args()
    pocket_number = args.pckt
    PDB_ID = args.pdb_id
    receptor_name = PDB_ID
    ligands_file = args.ligands
    tol = args.tol

    # Set up paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    receptor_folder = os.path.join(script_dir, receptor_name)
    if not os.path.exists(receptor_folder):
        os.makedirs(receptor_folder, exist_ok=True)

    # Set up logging
    log_file = os.path.join(receptor_folder, f"{receptor_name}_console_output.log")
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s %(message)s')

    logging.info(f"Tolerance (Å): {tol}")
    print(f"Docking tolerance set to {tol} Å.")

    # Folder for ligands
    ligands_folder = os.path.join(script_dir, 'ligands')

    # Output results file
    results_file = os.path.join(receptor_folder, f"{receptor_name}_results.txt")

    # Start processing
    try:
        # Download and prepare receptor
        downloaded_pdb_path = download_pdb(PDB_ID, receptor_folder)
        dirty_pdb = os.path.join(receptor_folder, f'{receptor_name}_dirty.pdb')
        shutil.move(downloaded_pdb_path, dirty_pdb)

        fixed_pdb = os.path.join(receptor_folder, f'{receptor_name}_fixed.pdb')
        fix_pdb(dirty_pdb, fixed_pdb, ph=7.4)

        receptor_pdbqt = os.path.join(receptor_folder, f"{receptor_name}.pdbqt")
        prepare_receptor(fixed_pdb, receptor_pdbqt)

        # Run P2Rank
        output_dir = os.path.join(receptor_folder, 'p2rank_output')
        run_p2rank(fixed_pdb, output_dir)

        # Get docking box parameters
        center_x, center_y, center_z, Size_x, Size_y, Size_z = get_docking_box(output_dir, fixed_pdb, tol, pocket_number)

        # Read ligands from SDF file
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

                # Prepare ligand
                prepare_ligand(ligand_pdb, ligand_pdbqt)

                # Run docking
                vina_output = run_vina(receptor_pdbqt, ligand_pdbqt, output_pdbqt,
                                    center_x, center_y, center_z,
                                    Size_x, Size_y, Size_z,
                                    exhaustiveness=16)

                # Save vina output to results file
                rf.write(f"Ligand: {ligand_name}\n")
                rf.write(vina_output)
                rf.write("\n\n")

                # Display information in the terminal after docking is complete
                print(f'Ligand {ligand_name} docked successfully!')

                # Generate visualization
                generate_visualizations(receptor_pdbqt, output_pdbqt, ligand_folder, receptor_name, ligand_name)

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        print(f"An error occurred: {e}")
        sys.exit(1)

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
    return expected_filename

@logger_decorator
def fix_pdb(input_pdb, output_pdb, ph=7.0):
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
        center_x = float(df['center_x'].iloc[0])
        center_y = float(df['center_y'].iloc[0])
        center_z = float(df['center_z'].iloc[0])

        # Selection of residuals for choosen pocket 
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
        return center_x, center_y, center_z, Size_x, Size_y, Size_z
    except Exception as e:
        logging.error(f"Error in getting docking box parameters: {e}")
        print(f"Error in getting docking box parameters: {e}")
        raise


@logger_decorator
def run_vina(receptor_pdbqt, ligand_pdbqt, output_pdbqt,
             center_x, center_y, center_z,
             size_x, size_y, size_z,
             exhaustiveness=8):
    try:
        vina_command = [
            VINA_PATH,
            '--receptor', receptor_pdbqt,
            '--ligand', ligand_pdbqt,
            '--out', output_pdbqt,
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

        return table_output

    except subprocess.CalledProcessError as e:
        logging.error(f"Error in docking with Vina: {e}")
        print(f"Error in docking with Vina: {e}")
        raise

@logger_decorator
def generate_visualizations(receptor_pdbqt, output_pdbqt, output_folder, receptor_name, ligand_name):
    try:
        cmd.load(output_pdbqt, 'docked_ligand')
        cmd.load(receptor_pdbqt, 'receptor')
        cmd.color('cyan', 'receptor')
        cmd.color('red', 'docked_ligand')
        set_visualization_and_focus()
        image_path = os.path.join(output_folder, f"{receptor_name}_{ligand_name}_docking.png")
        cmd.ray(1920, 1080)
        cmd.png(image_path, width=1920, height=1080, dpi=300)
        logging.info(f"Visualization saved: {image_path}")
        cmd.delete('all')
    except Exception as e:
        logging.error(f"Error in generating visualization: {e}")
        print(f"Error in generating visualization: {e}")
        raise

def set_visualization_and_focus():
    cmd.show('cartoon', 'all')
    cmd.bg_color('white')
    cmd.orient()
    cmd.zoom('all', buffer=1)
    cmd.clip('near', -10)
    cmd.clip('far', 10)
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

if __name__ == "__main__":
    main()