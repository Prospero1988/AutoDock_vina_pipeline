#!/usr/bin/env python3
"""
This script automates the process of docking multiple ligands to multiple receptor proteins using AutoDock Vina.
It reads receptors PDB ID from a CSV file, downloads each receptor structure from the PDB database, prepares them for docking, and processes multiple ligands from an SDF or MOL2 file.

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

- --pdb_ids: The name of the CSV file containing PDB IDs of receptor proteins, located in ./receptors directory.
- --ligands: The name of the SDF or MOL2 file containing ligands to dock, located in ./ligands directory.
- --tol: Tolerance in Ångströms to expand the docking pocket dimensions beyond those defined by P2Rank (default: 0).
- --pckt: Pocket number to use from P2Rank predictions (default: 1).
- --exhaust: Specifies how thorough the search should be for the best binding poses.
              Higher values increase precision but require more computation time (default: 16).
- --energy_range: Determines the range of energy scores (in kcal/mol) for poses to be considered (default: 3).

All receptor-related files will be saved in the ./PDB_ID directory. Each ligand's docking results will
be saved in ./PDB_ID/02_ligands_results/ligand_name_or_number.

Ensure that all required tools and libraries are installed and properly configured in your environment.
"""

import argparse
import os
import sys
import shutil
import subprocess
import logging
import re
import tempfile
import csv  # Import csv module

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
from pymol.cgo import CYLINDER
from pymol.cgo import BEGIN, LINES, VERTEX, END, COLOR, LINEWIDTH

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
    parser.add_argument('--ligands', required=True, help='Name of the SDF or MOL2 file containing ligands, located in ./ligands.')
    parser.add_argument('--tol_x', type=int, default=0, help='Tolerance in Ångströms to expand the docking pocket dimension in X beyond those defined by P2Rank (default: 0).')
    parser.add_argument('--tol_y', type=int, default=0, help='Tolerance in Ångströms to expand the docking pocket dimension in Y beyond those defined by P2Rank (default: 0).')
    parser.add_argument('--tol_z', type=int, default=0, help='Tolerance in Ångströms to expand the docking pocket dimension in Z beyond those defined by P2Rank (default: 0).')
    parser.add_argument('--pckt', type=int, default=1, help='Pocket number to use from P2Rank predictions (default: 1).')
    parser.add_argument('--exhaust', type=int, default=16, help='Specifies how thorough the search should be for the best binding poses. Higher values increase precision but require more computation time (default: 16).')
    parser.add_argument('--energy_range', type=int, default=3, help='Determines the range of energy scores (in kcal/mol) for poses to be considered (default: 3).')

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

    molecules = []

    if ligands_path.endswith(".sdf"):
        logging.info("Detected SDF file. Reading ligands...")
        suppl = Chem.SDMolSupplier(ligands_path)
        if not suppl:
            raise ValueError(f"Could not read ligands from {ligands_path} or file is empty.")
        molecules = [mol for mol in suppl if mol is not None]

    elif ligands_path.endswith(".mol2"):
        logging.info("Detected MOL2 file. Converting to SDF...")
        # Create a temporary SDF file
        with tempfile.NamedTemporaryFile(delete=False, suffix='.sdf') as tmp_sdf:
            tmp_sdf_path = tmp_sdf.name

        try:
            # Convert MOL2 to SDF using Open Babel
            conversion_command = [
                OBABEL_PATH,
                "-i", "mol2",
                ligands_path,
                "-o", "sdf",
                "-O", tmp_sdf_path,
                "--aromatic"
            ]
            subprocess.run(conversion_command, check=True)
            logging.info(f"Converted MOL2 to SDF: {tmp_sdf_path}")

            # Read the SDF file
            suppl = Chem.SDMolSupplier(tmp_sdf_path)
            if not suppl:
                raise ValueError(f"Could not read ligands from {tmp_sdf_path} or file is empty.")
            molecules = [mol for mol in suppl if mol is not None]
        except subprocess.CalledProcessError as e:
            logging.error(f"Error in converting MOL2 to SDF: {e}")
            print(f"Error in converting MOL2 to SDF: {e}")
            raise
        finally:
            # Remove the temporary SDF file
            if os.path.exists(tmp_sdf_path):
                os.remove(tmp_sdf_path)
                logging.info(f"Temporary SDF file removed: {tmp_sdf_path}")

    else:
        logging.error(f"Unsupported file format: {ligands_path}. Supported formats are .sdf and .mol2.")
        raise ValueError(f"Unsupported file format: {ligands_path}. Please provide .sdf or .mol2 files.")

    if not molecules:
        logging.warning(f"No valid molecules found in {ligands_path}.")
        raise ValueError(f"No valid molecules found in {ligands_path}.")
    else:
        logging.info(f"Number of valid molecules read: {len(molecules)}")

    # Assign names to molecules with exception handling and extract properties
    ligand_data_list = []
    for idx, mol in enumerate(molecules):
        if mol.HasProp('_Name') and mol.GetProp('_Name').strip() and not all(c == '*' for c in mol.GetProp('_Name').strip()):
            ligand_name = sanitize_ligand_name(mol.GetProp('_Name').strip())
        else:
            ligand_name = f"ligand_{idx + 1:03d}"  # Format: ligand_001, ligand_002, etc.
            logging.info(f"Assigned default name to molecule {idx + 1}: {ligand_name}")
        mol.SetProp('_Name', ligand_name)

        # Extract 'Code name' property if it exists
        if mol.HasProp('Code name'):
            code_name = mol.GetProp('Code name').strip()
        else:
            code_name = None

        # Extract 'Product (SMILES)' property if it exists
        if mol.HasProp('Product (SMILES)'):
            smiles = mol.GetProp('Product (SMILES)').strip()
        else:
            # Generate SMILES code from mol
            smiles = Chem.MolToSmiles(mol)

        ligand_data_list.append({
            'mol': mol,
            'name': ligand_name,
            'code_name': code_name,
            'smiles': smiles
        })

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
        tol_x = args.tol_x
        tol_y = args.tol_y
        tol_z = args.tol_z
        exhaustiveness = args.exhaust
        energy_range = args.energy_range

        # Set up paths
        receptor_folder = os.path.join(script_dir, receptor_name)
        if not os.path.exists(receptor_folder):
            os.makedirs(receptor_folder, exist_ok=True)

        # Creation of additional folder 03_ligands_PDBQT
        ligands_pdbqt_folder = os.path.join(receptor_folder, '03_ligands_PDBQT')
        os.makedirs(ligands_pdbqt_folder, exist_ok=True)
        logging.info(f"Created folder for ligands PDBQT: {ligands_pdbqt_folder}")

        # Creation of additional folder 02_ligands_results
        ligands_results_folder = os.path.join(receptor_folder, '02_ligands_results')
        os.makedirs(ligands_results_folder, exist_ok=True)
        logging.info(f"Created folder for ligands results: {ligands_results_folder}")

        # Set up logging for this receptor
        log_file = os.path.join(receptor_folder, f"{receptor_name}_console_output.log")
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s %(message)s')
        file_handler.setFormatter(formatter)
        logging.getLogger().addHandler(file_handler)

        logging.info(f"Processing receptor {PDB_ID}")
        logging.info(f"Tolerances (Å): X: {tol_x}, Y: {tol_y}, Z: {tol_z}")
        print(f"Docking tolerances set to X: {tol_x} Å, Y: {tol_y} Å, Z: {tol_z} Å.")

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
            center_x, center_y, center_z, Size_x, Size_y, Size_z, predictions_csv = get_docking_box(output_dir, fixed_pdb, tol_x, tol_y, tol_z, pocket_number)

            with open(results_file, 'w') as rf:
                for ligand_data in ligand_data_list:
                    mol = ligand_data['mol']
                    ligand_name = ligand_data['name']
                    code_name = ligand_data['code_name']
                    smiles = ligand_data['smiles']

                    if mol is None:
                        logging.warning(f"Skipping invalid molecule: {ligand_name}")
                        continue

                    # Placing ligands folders in 02_ligands_results
                    ligand_folder = os.path.join(ligands_results_folder, ligand_name)
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

                    # Copy PDBQT file after docking to folder 03_ligands_PDBQT
                    shutil.copy2(output_pdbqt, ligands_pdbqt_folder)
                    logging.info(f"Copied {output_pdbqt} to {ligands_pdbqt_folder}")

                    # Save vina output to results file
                    rf.write(f"Ligand: {ligand_name}\n")
                    rf.write(vina_output)
                    rf.write("\n\n")

                    # Display information in the terminal after docking is complete
                    print(f'Ligand {ligand_name} docked successfully!')

                    # Generate visualization
                    generate_visualizations(
                        receptor_pdbqt, output_pdbqt, ligand_folder, receptor_name, ligand_name,
                        center=(center_x, center_y, center_z),
                        size=(Size_x, Size_y, Size_z)
                    )

                    docking_image_path = os.path.join(ligand_folder, f"{receptor_name}_{ligand_name}_docking.png")

                    # Get the first affinity value
                    if affinities:
                        affinity = affinities[0][0]
                    else:
                        affinity = None

                    pymol_session_path = os.path.join(ligand_folder, f"{receptor_name}_{ligand_name}_docking.pse")

                    # Collect ligand results
                    ligand_results.append({
                        'name': ligand_name,
                        'image': image_filename,
                        'affinity': affinity,
                        'output_pdbqt': output_pdbqt,
                        'docking_image': docking_image_path,
                        'smiles': smiles,
                        'code_name': code_name,
                        'pymol_session': pymol_session_path
                    })

            # Generate HTML results file
            html_file = os.path.join(receptor_folder, f"{receptor_name}_results.html")
            generate_html_results(html_file, receptor_name, ligands_file, ligand_results, predictions_csv, protein_name, args.pckt, receptor_pdbqt)

            # Generate CSV results file
            csv_file = os.path.join(receptor_folder, f"{receptor_name}_results_in_CSV.csv")
            try:
                # Determine if any ligand has a code_name
                include_code_name = any(ligand['code_name'] is not None for ligand in ligand_results)

                # Open CSV file for writing
                with open(csv_file, 'w', newline='', encoding='utf-8') as csvfile:
                    fieldnames = ['name', 'affinity', 'smiles']
                    if include_code_name:
                        fieldnames.append('code_name')

                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()

                    for ligand in ligand_results:
                        affinity = ligand['affinity']
                        if affinity is not None:
                            affinity_str = f"{affinity:.2f}"
                        else:
                            affinity_str = ''
                        row = {
                            'name': ligand['name'],
                            'affinity': affinity_str,
                            'smiles': ligand['smiles'] if ligand['smiles'] else ''
                        }
                        if include_code_name:
                            row['code_name'] = ligand['code_name'] if ligand['code_name'] else ''
                        writer.writerow(row)

                logging.info(f"CSV results file generated: {csv_file}")
            except Exception as e:
                logging.error(f"Error in generating CSV results: {e}")
                print(f"Error in generating CSV results: {e}")
                raise

            print(f'Results HTML file generated: {html_file}')

        except Exception as e:
            logging.error(f"An error occurred with PDB ID {PDB_ID}: {e}")
            print(f"An error occurred with PDB ID {PDB_ID}: {e}")
            # Remove the file handler after processing each receptor
            logging.getLogger().removeHandler(file_handler)
            continue  # Proceed to the next PDB_ID

        # Remove the file handler after processing each receptor
        logging.getLogger().removeHandler(file_handler)

# Function to sanitize ligand names
def sanitize_ligand_name(name):
    """
    Sanitize the ligand name by replacing invalid characters with underscores.
    """
    # Define a regex pattern for valid characters (alphanumeric and underscores)
    valid_pattern = re.compile(r'[^A-Za-z0-9_-]')
    sanitized_name = valid_pattern.sub('_', name)
    # Additionally, remove leading/trailing underscores
    sanitized_name = sanitized_name.strip('_')
    # If the sanitized name is empty, assign a default name
    if not sanitized_name:
        sanitized_name = "ligand_unnamed"
    return sanitized_name

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
def get_docking_box(output_dir, receptor_pdb, tol_x, tol_y, tol_z, pocket_number):
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
        Size_x, Size_y, Size_z = cube_size + np.array([tol_x, tol_y, tol_z])

        logging.info(f"Docking box dimensions before tolerance: {cube_size}")
        logging.info(f"Docking box dimensions after tolerance: {Size_x}, {Size_y}, {Size_z}")
        logging.info(f"Tolerances applied: X: {tol_x}, Y: {tol_y}, Z: {tol_z}")
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

def add_axes(center=(0, 0, 0), length=10.0, radius=0.1, offset=(15, 15, 15)):
    """
    Add XYZ axes at a specified center position in PyMOL with labels.

    Parameters:
        center (tuple): Coordinates (x, y, z) for the center of the axes.
        length (float): Length of each axis in Ångströms.
        radius (float): Radius of the axes cylinders.
        offset (tuple): Offset to apply to the center coordinates.
    """

    x, y, z = center
    ox, oy, oz = offset
    new_center = (x + ox, y + oy, z + oz)

    axes = [
        # X-axis (red)
        CYLINDER, new_center[0], new_center[1], new_center[2], new_center[0] + length, new_center[1], new_center[2], radius,
        1.0, 0.0, 0.0,  # Kolor początkowy (czerwony)
        1.0, 0.0, 0.0,  # Kolor końcowy (czerwony)
        # Y-axis (green)
        CYLINDER, new_center[0], new_center[1], new_center[2], new_center[0], new_center[1] + length, new_center[2], radius,
        0.0, 1.0, 0.0,  # Kolor początkowy (zielony)
        0.0, 1.0, 0.0,  # Kolor końcowy (zielony)
        # Z-axis (blue)
        CYLINDER, new_center[0], new_center[1], new_center[2], new_center[0], new_center[1], new_center[2] + length, radius,
        0.0, 0.0, 1.0,  # Kolor początkowy (niebieski)
        0.0, 0.0, 1.0   # Kolor końcowy (niebieski)
    ]

    # Załaduj osie jako obiekt CGO
    cmd.load_cgo(axes, "axes")

    # Tworzenie pseudoatomów na końcach osi
    cmd.pseudoatom(object='axis_x_end', pos=(new_center[0] + length, new_center[1], new_center[2]))
    cmd.pseudoatom(object='axis_y_end', pos=(new_center[0], new_center[1] + length, new_center[2]))
    cmd.pseudoatom(object='axis_z_end', pos=(new_center[0], new_center[1], new_center[2] + length))

    # Etykietowanie pseudoatomów
    cmd.label('axis_x_end', '"X"')
    cmd.label('axis_y_end', '"Y"')
    cmd.label('axis_z_end', '"Z"')

    # Opcjonalnie, ukrycie pseudoatomów
    # cmd.hide('everything', 'axis_x_end axis_y_end axis_z_end')


    # Optionally, hide the pseudoatoms (uncomment if desired)
    #cmd.hide('everything', 'axis_x_end axis_y_end axis_z_end')


@logger_decorator
def generate_visualizations(receptor_pdbqt, output_pdbqt, output_folder, receptor_name, ligand_name, center, size, offset=(15, 15, 15)):
    try:
        # Załaduj struktury do PyMOL z plików PDBQT (zawierających ładunki)
        cmd.load(receptor_pdbqt, 'receptor')
        cmd.load(output_pdbqt, 'ligand')

        # Dodaj atomy wodoru (jeśli nie są obecne)
        cmd.h_add('receptor')
        cmd.h_add('ligand')

        # Ustaw kolory receptor i ligandu
        cmd.color('cyan', 'receptor')
        cmd.color('red', 'ligand')

        # Ustaw reprezentacje
        cmd.show('cartoon', 'receptor')
        cmd.show('sticks', 'ligand')

        # Ukryj atomy wodoru w wizualizacji
        cmd.hide('sticks', '(elem H)')

        # Dodaj osie współrzędnych z przesunięciem
        add_axes(center=center, length=10.0, radius=0.3, offset=offset)

        # Rysuj box dokowania
        draw_docking_box(center, size)

        # Dostosuj wizualizację i ustaw kamerę
        set_visualization_and_focus()

        # Zapisz obrazek jako PNG
        image_path = os.path.join(output_folder, f"{receptor_name}_{ligand_name}_docking.png")
        cmd.ray(1920, 1080)
        cmd.png(image_path, width=1920, height=1080, dpi=300)
        logging.info(f"Visualization saved: {image_path}")

        # Zapisz sesję PyMOL
        session_path = os.path.join(output_folder, f"{receptor_name}_{ligand_name}_docking.pse")
        cmd.save(session_path)
        logging.info(f"PyMOL session saved: {session_path}")

        # Wyczyść sesję PyMOL
        cmd.delete('all')
    except Exception as e:
        logging.error(f"Error in generating visualization: {e}")
        print(f"Error in generating visualization: {e}")
        raise


def draw_docking_box(center, size):
    x, y, z = center
    sx, sy, sz = size
    # Calculate the min and max coordinates
    min_x = x - sx / 2
    max_x = x + sx / 2
    min_y = y - sy / 2
    max_y = y + sy / 2
    min_z = z - sz / 2
    max_z = z + sz / 2

    # Create the box using CGO objects
    
    box = [
        LINEWIDTH, 1.0,
        COLOR, 1.0, 1.0, 0.0,  # Yellow color
        BEGIN, LINES,
        VERTEX, min_x, min_y, min_z,
        VERTEX, max_x, min_y, min_z,

        VERTEX, max_x, min_y, min_z,
        VERTEX, max_x, max_y, min_z,

        VERTEX, max_x, max_y, min_z,
        VERTEX, min_x, max_y, min_z,

        VERTEX, min_x, max_y, min_z,
        VERTEX, min_x, min_y, min_z,

        VERTEX, min_x, min_y, max_z,
        VERTEX, max_x, min_y, max_z,

        VERTEX, max_x, min_y, max_z,
        VERTEX, max_x, max_y, max_z,

        VERTEX, max_x, max_y, max_z,
        VERTEX, min_x, max_y, max_z,

        VERTEX, min_x, max_y, max_z,
        VERTEX, min_x, min_y, max_z,

        VERTEX, min_x, min_y, min_z,
        VERTEX, min_x, min_y, max_z,

        VERTEX, max_x, min_y, min_z,
        VERTEX, max_x, min_y, max_z,

        VERTEX, max_x, max_y, min_z,
        VERTEX, max_x, max_y, max_z,

        VERTEX, min_x, max_y, min_z,
        VERTEX, min_x, max_y, max_z,
        END
    ]

    cmd.load_cgo(box, 'docking_box')


def set_visualization_and_focus():
    cmd.show('cartoon', 'all')
    cmd.bg_color('white')
    cmd.center('ligand')  # Centrowanie na ligandzie
    cmd.zoom('ligand', buffer=1.0)  # Powiększenie na ligandę
    cmd.move('z', -50)  # Przesunięcie kamery wzdłuż osi Z

    # Usunięcie ustawień przycinania
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

@logger_decorator
def generate_html_results(html_file, receptor_name, ligands_file, ligand_results, predictions_csv, protein_name, pckt, receptor_pdbqt):
    try:
        # Wczytaj dane z P2Rank i usuń nadmiarowe spacje
        p2rank_csv = predictions_csv
        df_p2rank = pd.read_csv(p2rank_csv)
        df_p2rank.columns = df_p2rank.columns.str.strip()
        df_p2rank = df_p2rank.map(lambda x: x.strip() if isinstance(x, str) else x)
        df_p2rank['score'] = df_p2rank['score'].astype(float).map("{:.2f}".format)
        df_p2rank['probability'] = df_p2rank['probability'].astype(float).map("{:.2f}".format)

        # Sortowanie wyników ligandów według 'affinity' od najniższej do najwyższej
        # Ligandy z affinity=None zostaną umieszczone na końcu
        ligand_results_sorted = sorted(
            ligand_results,
            key=lambda x: (x['affinity'] is None, x['affinity'] if x['affinity'] is not None else float('inf'))
        )

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
            hf.write('.docking-energy { background-color: #dfffe0; } /* Pastel green */\n')  # Stylizacja dla energii dokowania
            # Stylizacja dla kolumny "Docking Results" w pierwszej tabeli
            hf.write('td:nth-child(6), th:nth-child(6) {\n')
            hf.write('  max-width: 200px;\n')  # Maksymalna szerokość
            hf.write('  word-wrap: break-word;\n')  # Zawijanie tekstu
            hf.write('  white-space: normal;\n')  # Normalna biała spacja
            hf.write('  text-align: left;\n')  # Wyrównanie do lewej
            hf.write('}\n')
            # Stylizacja dla kolumny "score" w drugiej tabeli z dodatkowym marginesem
            hf.write('.p2rank-table td:nth-child(3), .p2rank-table th:nth-child(3) {\n')
            hf.write('  white-space: nowrap;\n')  # Dopasowanie szerokości do zawartości
            hf.write('  padding-left: 15px;\n')  # Lewy padding
            hf.write('  padding-right: 15px;\n')  # Prawy padding
            hf.write('}\n')
            # Stylizacja dla kolumny "residue_ids" w drugiej tabeli
            hf.write('.p2rank-table td:nth-child(5), .p2rank-table th:nth-child(5) {\n')
            hf.write('  max-width: 400px;\n')  # Ograniczenie szerokości kolumny
            hf.write('  word-wrap: break-word;\n')  # Zawijanie tekstu w kolumnie
            hf.write('  white-space: normal;\n')  # Normalna biała spacja dla zawijania
            hf.write('  padding: 15px;\n')  # Wewnętrzne marginesy
            hf.write('  text-align: left;\n')  # Opcjonalnie: wyrównanie tekstu do lewej
            hf.write('}\n')
            hf.write('</style>\n')

            hf.write('</head>\n')
            hf.write('<body>\n')

            # Nagłówek tabeli z wynikami dokowania
            receptor_pdbqt_name = os.path.basename(receptor_pdbqt)
            header_text = (
                f'Docking results for receptor with PDB code: '
                f'<span style="color: red;">{receptor_name}</span></br>'
                f'using structures from file: '
                f'<span style="color: navy;">{ligands_file}</span></br>'
            )
            hf.write(f'<h2 style="text-align: center;">{header_text}</h2>\n')
            hf.write(f'</br><h2 style="text-align: center; color: green; max-width: 600px; word-wrap: break-word; margin: auto;">{protein_name}</h2></br>')
            hf.write(f'<h3 style="text-align: center;">Docking to pocket number: <span style="color: red;">{pckt}</span></br></h2>\n')
            hf.write(f'<div style="text-align: center;">')
            hf.write(
                f'<a href="{receptor_pdbqt}" download="receptor_structure.pdbqt" type="application/octet-stream">'
                f'Receptor structure (file .PDBQT). DOWNLOAD</a>'
            )
            hf.write('</div>')
            hf.write('</br>')

            # Pierwsza tabela: Wyniki dokowania (posortowane według energii dokowania)
            hf.write('<table>\n')
            hf.write(
                '<tr>'
                '<th>Number</th>'
                '<th>Compound Name</th>'
                '<th>Structure</th>'
                '<th>Docking Image</th>'
                '<th class="docking-energy">Docking Energy<br/>(kcal/mol)</th>'
                '<th>Docking Results</th>'
                '<th>PyMOL Session</th>'
                '</tr>\n'
            )
            for idx, result in enumerate(ligand_results_sorted, start=1):
                name = result['name']
                image_path = os.path.relpath(result['image'], os.path.dirname(html_file))
                docking_image_path = os.path.relpath(result['docking_image'], os.path.dirname(html_file))
                affinity = result['affinity']
                if affinity is not None:
                    affinity_str = f"{affinity:.2f}"
                else:
                    affinity_str = 'N/A'
                # Ścieżka do pliku wynikowego PDBQT
                output_pdbqt_path = os.path.relpath(result['output_pdbqt'], os.path.dirname(html_file))
                link_text = os.path.basename(result['output_pdbqt'])
                # Ścieżka do pliku sesji PyMOL
                pymol_session_path = os.path.relpath(result['pymol_session'], os.path.dirname(html_file))
                pymol_session_link = os.path.basename(result['pymol_session'])

                hf.write('<tr>\n')
                hf.write(f'<td>{idx}</td>\n')  # 1. kolumna: Numer
                hf.write(f'<td>{name}</td>\n')  # 2. kolumna: Nazwa związku
                hf.write(f'<td><img src="{image_path}" alt="{name}" width="400"/></td>\n')  # 3. kolumna: Struktura (SVG)
                # 4. kolumna: Obrazek z dokowaniem
                hf.write(
                    f'<td><a href="{docking_image_path}" target="_blank">'
                    f'<img src="{docking_image_path}" alt="Docking Image" style="max-width: 150px; max-height: 150px;"></a></td>\n'
                )
                hf.write(f'<td class="docking-energy">{affinity_str}</td>\n')  # 5. kolumna: Energia dokowania
                # 6. kolumna: Wyniki dokowania (plik PDBQT)
                hf.write(
                    f'<td><a href="{output_pdbqt_path}" download="{link_text}" type="application/octet-stream">'
                    f'{link_text}</a></td>\n'
                )
                # 7. kolumna: Sesja PyMOL (.pse)
                hf.write(
                    f'<td><a href="{pymol_session_path}" download="{pymol_session_link}" '
                    f'type="application/octet-stream">{pymol_session_link}</a></td>\n'
                )
                hf.write('</tr>\n')
            hf.write('</table>\n')
            hf.write('</br>')

            # Link do szczegółowych wyników
            results_file = f"{receptor_name}_results.txt"
            hf.write('<p style="text-align: center; margin-top: 20px;">\n')
            hf.write(
                f'<a href="{results_file}" target="_blank">'
                f'Detailed results for each compound (All docking poses). CLICK</a>\n'
            )
            hf.write('</p>\n')
            hf.write('</br></br>')

            # Nagłówek dla tabeli z danymi P2Rank
            p2rank_header = f'P2RANK: identified docking pockets for receptor with PDB code: {receptor_name}'
            hf.write(f'<h3 style="text-align: center; margin-top: 20px;">{p2rank_header}</h3>\n')

            # Druga tabela: Dane z P2Rank
            hf.write('<table class="p2rank-table">\n')  # Dodanie klasy "p2rank-table"
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
            # Link do szczegółowych wyników P2Rank
            residues_csv_file = os.path.join('01_p2rank_output', f'{receptor_name}_fixed.pdb_residues.csv')
            hf.write('<p style="text-align: center; margin-top: 20px;">\n')
            hf.write(
                f'<a href="{residues_csv_file}" target="_blank">'
                f'Detailed information about individual amino acids </br>and their involvement in docking pockets. CLICK</a>\n'
            )
            hf.write('</p>\n')

            hf.write('</br></br>')

            # Informacje o autorze na końcu
            hf.write('<p style="font-size: small; text-align: center; margin-top: 20px;">\n')
            hf.write('Docking system based on <b>AutoDock Vina v.1.2.5</b> and <b>P2RANK v.2.4.2</b><br/>\n')
            hf.write('<b>Author:</b> Arkadiusz Leniak <b>email:</b> arkadiusz.leniak@gmail.com<br/>\n')
            hf.write(
                '<b>github:</b> '
                '<a href="https://github.com/Prospero1988/AutoDock_vina_pipeline">'
                'https://github.com/Prospero1988/AutoDock_vina_pipeline</a>\n'
            )
            hf.write('</p>\n')

            hf.write('</body>\n')
            hf.write('</html>\n')
        logging.info(f"HTML results file generated: {html_file}")
    except Exception as e:
        logging.error(f"Error in generating HTML results: {e}")
        print(f"Error in generating HTML results: {e}")
        raise


def sanitize_ligand_name(name):
    """
    Sanitize the ligand name by replacing invalid characters with underscores.
    """
    # Define a regex pattern for valid characters (alphanumeric and underscores)
    valid_pattern = re.compile(r'[^A-Za-z0-9_-]')
    sanitized_name = valid_pattern.sub('_', name)
    # Additionally, remove leading/trailing underscores
    sanitized_name = sanitized_name.strip('_')
    # If the sanitized name is empty, assign a default name
    if not sanitized_name:
        sanitized_name = "ligand_unnamed"
    return sanitized_name

if __name__ == "__main__":
    main()
