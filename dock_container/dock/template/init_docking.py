#!/usr/bin/env python3

import argparse
import os
import sys
import shutil
import subprocess
import logging
import re
import tempfile
import csv
from pathlib import Path
import gzip

from Bio import PDB
from Bio.PDB import PDBList
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
import pandas as pd
import numpy as np
import pymol
from pymol import cmd
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
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

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Automated docking script for multiple ligands and receptors.")
    parser.add_argument('--pdb_ids', required=True, help='Name of the CSV file containing PDB IDs and Chain IDs of receptor proteins, located in ./receptors.')
    parser.add_argument('--ligands', required=True, help='Name of the SDF or MOL2 file containing ligands, located in ./ligands.')
    parser.add_argument('--tol_x', type=int, default=0, help='Tolerance in Å to expand the docking pocket dimension in X (default: 0).')
    parser.add_argument('--tol_y', type=int, default=0, help='Tolerance in Å to expand the docking pocket dimension in Y (default: 0).')
    parser.add_argument('--tol_z', type=int, default=0, help='Tolerance in Å to expand the docking pocket dimension in Z (default: 0).')
    parser.add_argument('--offset_x', type=float, default=0.0, help='Offset in Å to shift the center of the docking grid box along X-axis (default: 0).')
    parser.add_argument('--offset_y', type=float, default=0.0, help='Offset in Å along Y-axis (default: 0).')
    parser.add_argument('--offset_z', type=float, default=0.0, help='Offset in Å along Z-axis (default: 0).')
    parser.add_argument('--pckt', type=int, default=1, help='Pocket number to use from P2Rank predictions (default: 1).')
    parser.add_argument('--exhaust', type=int, default=16, help='Exhaustiveness of the search (default: 16).')
    parser.add_argument('--energy_range', type=int, default=3, help='Energy range in kcal/mol (default: 3).')
    parser.add_argument('--num_modes', type=int, default=20, help='Number of maximum generated conformers (default: 20).')
    parser.add_argument('--seed', type=int, default=1988, help='Seed for the random number generator (default: 1988).')
    parser.add_argument('--flex', default=None, help='Path to flexible receptor PDBQT file (optional).')
    parser.add_argument('--rigid', default=None, help='Path to rigid receptor PDBQT file (optional).')
    parser.add_argument('--relax', action='store_true', help='Perform structure relaxation with MM on UFF with only 5 steps (optional).')
    parser.add_argument('--keepids', action='store_true', help='Retain original chain ID and residue numbering.')

    args = parser.parse_args()

    if args.flex == 'None':
        args.flex = None
    if args.rigid == 'None':
        args.rigid = None

    seed = args.seed

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
            suppl = Chem.SDMolSupplier(tmp_sdf_path)
            if not suppl:
                raise ValueError(f"Could not read ligands from {tmp_sdf_path} or file is empty.")
            molecules = [mol for mol in suppl if mol is not None]
        except subprocess.CalledProcessError as e:
            logging.error(f"Error in converting MOL2 to SDF: {e}")
            print(f"Error in converting MOL2 to SDF: {e}")
            raise
        finally:
            if os.path.exists(tmp_sdf_path):
                os.remove(tmp_sdf_path)
                logging.info(f"Temporary SDF file removed: {tmp_sdf_path}")
    else:
        logging.error(f"Unsupported file format: {ligands_path}. Supported: .sdf, .mol2.")
        raise ValueError(f"Unsupported file format: {ligands_path}. Provide .sdf or .mol2 files.")

    if not molecules:
        logging.warning(f"No valid molecules found in {ligands_path}.")
        raise ValueError(f"No valid molecules found in {ligands_path}.")
    else:
        logging.info(f"Number of valid molecules read: {len(molecules)}")

    # Assign names to molecules and extract properties
    ligand_data_list = []
    for idx, mol in enumerate(molecules):
        if mol.HasProp('_Name') and mol.GetProp('_Name').strip() and not all(c == '*' for c in mol.GetProp('_Name').strip()):
            ligand_name = sanitize_ligand_name(mol.GetProp('_Name').strip())
        else:
            ligand_name = f"ligand_{idx + 1:03d}"
            logging.info(f"Assigned default name to molecule {idx + 1}: {ligand_name}")
        mol.SetProp('_Name', ligand_name)

        if mol.HasProp('Code name'):
            code_name = mol.GetProp('Code name').strip()
        else:
            code_name = None

        if mol.HasProp('Product (SMILES)'):
            smiles = mol.GetProp('Product (SMILES)').strip()
        else:
            smiles = Chem.MolToSmiles(mol)

        #mol = Chem.AddHs(mol)  # Dodaj brakujące atomy H
        #Chem.SanitizeMol(mol)
        
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

    df_receptors = pd.read_csv(receptors_file_path, header=None)
    pdb_ids_and_chains = df_receptors.values.tolist()

    if df_receptors.empty:
        raise ValueError("No PDB IDs and Chain IDs found in the receptors file.")

    for row in pdb_ids_and_chains:
        if len(row) < 2:
            logging.warning(f"Skipping row with insufficient columns: {row}")
            continue
        PDB_ID, chain_ID = row[:2]
        pocket_number = args.pckt
        receptor_name = f"{PDB_ID}_{chain_ID}"
        tol_x = args.tol_x
        tol_y = args.tol_y
        tol_z = args.tol_z
        exhaustiveness = args.exhaust
        energy_range = args.energy_range

        receptor_folder = os.path.join(script_dir, receptor_name)
        if not os.path.exists(receptor_folder):
            os.makedirs(receptor_folder, exist_ok=True)

        ligands_pdbqt_folder = os.path.join(receptor_folder, '03_ligands_PDBQT')
        os.makedirs(ligands_pdbqt_folder, exist_ok=True)
        logging.info(f"Created folder for ligands PDBQT: {ligands_pdbqt_folder}")

        ligands_results_folder = os.path.join(receptor_folder, '02_ligands_results')
        os.makedirs(ligands_results_folder, exist_ok=True)
        logging.info(f"Created folder for ligands results: {ligands_results_folder}")

        log_file = os.path.join(receptor_folder, f"{receptor_name}_console_output.log")
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s %(message)s')
        file_handler.setFormatter(formatter)
        logging.getLogger().addHandler(file_handler)

        logging.info(f"Processing receptor {PDB_ID} with Chain ID {chain_ID}")
        print(f"Processing receptor {PDB_ID} with Chain ID {chain_ID}")
        print(f"Docking tolerances set to X: {tol_x} Å, Y: {tol_y} Å, Z: {tol_z} Å.")

        results_file = os.path.join(receptor_folder, f"{receptor_name}_results.txt")
        ligand_results = []

        try:
            # Always download and prepare receptor
            downloaded_pdb_path, protein_name = download_pdb(PDB_ID, receptor_folder)
            dirty_pdb = os.path.join(receptor_folder, f'{receptor_name}_dirty.pdb')
            shutil.move(downloaded_pdb_path, dirty_pdb)

            fixed_pdb = os.path.join(receptor_folder, f'{receptor_name}_fixed.pdb')
            fix_pdb(dirty_pdb, fixed_pdb, chain_ID, args.relax, args.keepids, ph=7.4)

            receptor_pdbqt = os.path.join(receptor_folder, f"{receptor_name}.pdbqt")
            prepare_receptor(fixed_pdb, receptor_pdbqt)
            flex_pdbqt = None

            # Run P2Rank and get docking box
            output_dir = os.path.join(receptor_folder, '01_p2rank_output')
            run_p2rank(fixed_pdb, output_dir)
            center_x, center_y, center_z, Size_x, Size_y, Size_z, predictions_csv = get_docking_box(
                output_dir, fixed_pdb, tol_x, tol_y, tol_z, pocket_number,
                args.offset_x, args.offset_y, args.offset_z
            )

            # If user provided rigid and flex, override receptor_pdbqt and flex_pdbqt now
            if args.rigid is not None and args.flex is not None and args.rigid != 'None' and args.flex != 'None':
                logging.info("Flexible docking mode: using provided rigid and flex files.")
                receptor_pdbqt = os.path.abspath(args.rigid)
                flex_pdbqt = os.path.abspath(args.flex)

            with open(results_file, 'w') as rf:
                for ligand_data in ligand_data_list:
                    mol = ligand_data['mol']
                    ligand_name = ligand_data['name']
                    code_name = ligand_data['code_name']
                    smiles = ligand_data['smiles']

                    if mol is None:
                        logging.warning(f"Skipping invalid molecule: {ligand_name}")
                        continue

                    ligand_folder = os.path.join(ligands_results_folder, ligand_name)
                    os.makedirs(ligand_folder, exist_ok=True)

                    ligand_pdb = os.path.join(ligand_folder, f"{ligand_name}.pdb")
                    ligand_pdbqt = os.path.join(ligand_folder, f"{ligand_name}.pdbqt")
                    output_pdbqt = os.path.join(ligand_folder, f"{receptor_name}_{ligand_name}.pdbqt")

                    write_mol_to_pdb(mol, ligand_pdb)
                    AllChem.Compute2DCoords(mol)
                    image_filename = os.path.join(ligand_folder, f"{ligand_name}.svg")
                    draw_molecule_to_file(mol, image_filename)

                    prepare_ligand(ligand_pdb, ligand_pdbqt)

                    num_modes = args.num_modes
                    vina_output, affinities = run_vina(
                        receptor_pdbqt, ligand_pdbqt, output_pdbqt,
                        center_x, center_y, center_z,
                        Size_x, Size_y, Size_z,
                        exhaustiveness, energy_range, num_modes, seed, flex_pdbqt=flex_pdbqt
                    )

                    shutil.copy2(output_pdbqt, ligands_pdbqt_folder)
                    logging.info(f"Copied {output_pdbqt} to {ligands_pdbqt_folder}")

                    rf.write(f"Ligand: {ligand_name}\n")
                    rf.write(vina_output)
                    rf.write("\n\n")

                    print(f'Ligand {ligand_name} docked successfully!')

                    generate_visualizations(
                        receptor_pdbqt, output_pdbqt, ligand_folder, receptor_name, ligand_name,
                        center=(center_x, center_y, center_z),
                        size=(Size_x, Size_y, Size_z)
                    )

                    docking_image_path = os.path.join(ligand_folder, f"{receptor_name}_{ligand_name}_docking.png")

                    if affinities:
                        affinity = affinities[0][0]
                    else:
                        affinity = None

                    pymol_session_path = os.path.join(ligand_folder, f"{receptor_name}_{ligand_name}_docking.pse")

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

            html_file = os.path.join(receptor_folder, f"{receptor_name}_results.html")
            generate_html_results(html_file, receptor_name, ligands_file, ligand_results, predictions_csv, protein_name, args.pckt, receptor_pdbqt)

            csv_file = os.path.join(receptor_folder, f"{receptor_name}_results_in_CSV.csv")
            try:
                include_code_name = any(ligand['code_name'] is not None for ligand in ligand_results)
                with open(csv_file, 'w', newline='', encoding='utf-8') as csvfile:
                    fieldnames = ['name', 'affinity', 'smiles']
                    if include_code_name:
                        fieldnames.append('code_name')
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                    for ligand in ligand_results:
                        affinity = ligand['affinity']
                        affinity_str = f"{affinity:.2f}" if affinity is not None else ''
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
            logging.getLogger().removeHandler(file_handler)
            continue

        logging.getLogger().removeHandler(file_handler)

def sanitize_ligand_name(name):
    valid_pattern = re.compile(r'[^A-Za-z0-9_-]')
    sanitized_name = valid_pattern.sub('_', name)
    sanitized_name = sanitized_name.strip('_')
    if not sanitized_name:
        sanitized_name = "ligand_unnamed"
    return sanitized_name

@logger_decorator
def download_pdb(pdb_id, download_dir):
    """
    Próbuje ściągnąć plik PDB z serwera RCSB. Jeśli to się nie uda,
    pobiera plik mmCIF i konwertuje go do PDB za pomocą BioPython.
    """

    if not os.path.exists(download_dir):
        os.makedirs(download_dir)

    pdbl = PDBList()
    expected_pdb_filename = os.path.join(download_dir, f"{pdb_id}.pdb")

    try:
        logging.info(f"Attempting to download PDB file for {pdb_id}...")
        pdb_file_path = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb', pdir=download_dir)

        # Obsługa pliku .gz - wypakowanie
        if pdb_file_path.endswith('.gz'):
            with gzip.open(pdb_file_path, 'rb') as f_in:
                uncompressed_path = pdb_file_path[:-3]  # usuwamy .gz
                with open(uncompressed_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(pdb_file_path)
            pdb_file_path = uncompressed_path

        # Przenosimy na docelową nazwę (np. 1abc.pdb)
        shutil.move(pdb_file_path, expected_pdb_filename)
        logging.info(f"PDB file downloaded: {expected_pdb_filename}")

    except Exception as e_pdb:
        logging.warning(f"Failed to download PDB file for {pdb_id}: {e_pdb}")
        logging.info(f"Attempting to download mmCIF file for {pdb_id}...")

        try:
            cif_file_path = pdbl.retrieve_pdb_file(pdb_id, file_format='mmCif', pdir=download_dir)
            # Obsługa pliku .cif.gz
            if cif_file_path.endswith('.cif.gz'):
                with gzip.open(cif_file_path, 'rb') as f_in:
                    uncompressed_path = cif_file_path[:-3]
                    with open(uncompressed_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(cif_file_path)
                cif_file_path = uncompressed_path

            logging.info(f"mmCIF file downloaded: {cif_file_path}")
            logging.info(f"Converting mmCIF to PDB for {pdb_id} via BioPython...")

            # ---- OTO NAJWAŻNIEJSZA ZMIANA: KONWERSJA PRZEZ BioPython ----
            mmcif_parser = MMCIFParser(QUIET=True)
            structure = mmcif_parser.get_structure(pdb_id, cif_file_path)

            pdbio = PDBIO()
            pdbio.set_structure(structure)
            pdbio.save(expected_pdb_filename)
            logging.info(f"PDB file after conversion: {expected_pdb_filename}")

            # Sprzątamy plik .cif, jeśli już niepotrzebny
            if os.path.exists(cif_file_path):
                os.remove(cif_file_path)
                logging.info(f"Removed mmCIF file: {cif_file_path}")

        except Exception as e_cif:
            logging.error(f"Failed to download both PDB and mmCIF files for {pdb_id}: {e_cif}")
            raise FileNotFoundError(f"Could not download PDB or mmCIF file for {pdb_id}.")

    # Ostateczny check, czy mamy plik PDB:
    if not os.path.exists(expected_pdb_filename):
        raise FileNotFoundError(f"PDB file for {pdb_id} was not found after attempts.")

    # Opcjonalnie: wyciągamy nazwę białka z nagłówka PDB (może się nie udać)
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, expected_pdb_filename)
        header = structure.header
        protein_name = header.get('name', 'Unknown')
        logging.info(f"Protein name for {pdb_id}: {protein_name}")
    except Exception as e_header:
        logging.warning(f"Could not extract protein name for {pdb_id}: {e_header}")
        protein_name = 'Unknown'

    # Jeśli mamy protein_name, wstawiamy go w TITLE
    if protein_name != 'Unknown':
        with open(expected_pdb_filename, 'r') as f:
            pdb_lines = f.readlines()
        insert_index = 0
        for i, line in enumerate(pdb_lines):
            if line.startswith("HEADER"):
                insert_index = i + 1
                break
        title_line = f"TITLE     {protein_name}\n"
        pdb_lines.insert(insert_index, title_line)
        with open(expected_pdb_filename, 'w') as f:
            f.writelines(pdb_lines)
        logging.info("Inserted TITLE line with protein name into the PDB file.")

    return expected_pdb_filename, protein_name


@logger_decorator
def fix_pdb(input_pdb, output_pdb, chain_ID, relax, keepIds=False, ph=7.4):
    
    def pymol_fix_structure(input_pdb):
        """Naprawia strukturę PDB za pomocą PyMOL."""
        try:
            cmd.load(input_pdb, "molecule")
            cmd.remove("hydro")  # Usunięcie wszystkich wodorów
            cmd.h_add("polymer")  # Dodanie brakujących wodorów
            cmd.save(input_pdb)
            cmd.delete("all")
            logging.info(f"Structure fixed with PyMOL and saved to: {input_pdb}")
        except Exception as e:
            logging.error(f"Error during PyMOL structure fixing: {e}")
            raise
    
    
    try:
        fixer = PDBFixer(filename=input_pdb)
    
        # Downloading the list of strings
        chain_list = list(fixer.topology.chains())
        chain_ids = [chain.id for chain in chain_list]

        if chain_ID not in chain_ids:
            raise ValueError(f"Chain ID {chain_ID} not found in PDB structure {input_pdb}.")

        # Finding the index of the selected string
        chain_to_keep_index = chain_ids.index(chain_ID)

        logging.info(f"Maintaining a chain with an index {chain_to_keep_index} (ID '{chain_ID}')")
        print(f"Keeping the chain with the ID'{chain_ID}'")

        # Removal of all other chains
        chains_to_remove = [i for i in range(len(chain_list)) if i != chain_to_keep_index]
        fixer.removeChains(chains_to_remove)

        # Removal of heterogeneities (including water)
        fixer.removeHeterogens(keepWater=False)

        # Finding the missing residuals
        fixer.findMissingResidues()
        fixer.missingResidues = {}
        # Finding the missing atoms and adding them
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        # Adding missing protons
        fixer.addMissingHydrogens(ph)

        # Saving the repaired PDB file
        with open(output_pdb, 'w') as outfile:
            PDBFile.writeFile(fixer.topology, fixer.positions, outfile, keepIds=keepIds)
        logging.info(f"Fixed PDB saved as {output_pdb}")

        if relax:
            try:
                logging.debug("Fixing structure with PyMOL.")
                pymol_fix_structure(output_pdb)
                logging.debug("Relaxing structure with RDKit.")
                mol = Chem.MolFromPDBFile(output_pdb, removeHs=False)

                if mol is not None:
                    
                    # Embedding konformerów, jeśli brak
                    if mol.GetNumConformers() == 0:
                        try:
                            logging.debug("No conformers found, embedding molecule.")
                            AllChem.EmbedMolecule(mol)
                        except Exception as embed_error:
                            logging.warning(f"Error during conformer embedding: {embed_error}. Proceeding without embedding.")

                    # Próba optymalizacji UFF
                    try:
                        logging.debug("Performing UFF optimization (5 steps).")
                        AllChem.UFFOptimizeMolecule(mol, maxIters=5)
                        logging.info("UFF optimization completed.")
                    except Exception as uff_error:
                        logging.warning(f"Error during UFF optimization: {uff_error}. Trying MMFF.")

                        # Próba optymalizacji MMFF jako alternatywa
                        try:
                            logging.debug("Attempting MMFF optimization.")
                            if Chem.MMFFHasAllMoleculeParams(mol):
                                mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
                                AllChem.MMFFOptimizeMolecule(mol, mmff_props, maxIters=5)
                                logging.info("MMFF optimization completed.")
                            else:
                                logging.warning("MMFF parameters are not available for this molecule. Skipping MMFF optimization.")
                        except Exception as mmff_error:
                            logging.warning(f"Error during MMFF optimization: {mmff_error}. Proceeding without optimization.")

                    # Zapisanie zrelaksowanej struktury do pliku PDB
                    try:
                        logging.debug("Saving relaxed structure back to PDB.")
                        Chem.MolToPDBFile(mol, output_pdb)
                        logging.info("Relaxed structure saved successfully.")
                    except Exception as save_error:
                        logging.warning(f"Error saving relaxed PDB: {save_error}.")
                else:
                    logging.warning("RDKit could not read the fixed PDB for optimization.")
            except Exception as e:
                logging.error(f"Unexpected error in RDKit processing: {e}. Proceeding without RDKit relaxation.")   
    

        # Checking whether the file contains atoms
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('fixed', output_pdb)
        atom_count = len(list(structure.get_atoms()))
        if atom_count == 0:
            raise ValueError("Fixed PDB has no atoms after processing.")
        else:
            logging.info(f"Fixed PDB contains {atom_count} atoms.")
            print(f"Fixed PDB contains {atom_count} atoms.")

    except Exception as e:
        logging.error(f"Error in fixing PDB: {e}")
        print(f"Error in fixing PDB: {e}")
        raise

@logger_decorator
def prepare_receptor(input_pdb, output_pdbqt):
    try:
        subprocess.run(
            [OBABEL_PATH, "-i", "pdb", input_pdb, "-o", "pdbqt", "-O", output_pdbqt, "-xr", "--partialcharge", "gasteiger"],
            check=True
        )
        logging.info(f"Receptor prepared: {output_pdbqt}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in preparing receptor: {e}")
        print(f"Error in preparing receptor: {e}")
        raise

@logger_decorator
def prepare_ligand(input_pdb, output_pdbqt):
    try:
        subprocess.run(
            [OBABEL_PATH, "-i", "pdb", input_pdb, "-o", "pdbqt", "-O", output_pdbqt, "-h", "--gen3d", "--partialcharge", "gasteiger"],
            check=True
        )
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
def get_docking_box(output_dir, receptor_pdb, tol_x, tol_y, tol_z, pocket_number, offset_x, offset_y, offset_z):
    try:
        receptor_base_name = os.path.basename(receptor_pdb).split('.')[0]
        predictions_csv = os.path.join(output_dir, f'{receptor_base_name}.pdb_predictions.csv')
        residues_csv = os.path.join(output_dir, f'{receptor_base_name}.pdb_residues.csv')

        df = pd.read_csv(predictions_csv)
        pred = pd.read_csv(residues_csv)
        df.columns = df.columns.str.strip()
        pred.columns = pred.columns.str.strip()

        if pocket_number <= 0 or pocket_number > len(df):
            raise ValueError(f"Pocket number {pocket_number} is out of range. Available: 1-{len(df)}.")

        center_x = float(df['center_x'].iloc[pocket_number - 1])
        center_y = float(df['center_y'].iloc[pocket_number - 1])
        center_z = float(df['center_z'].iloc[pocket_number - 1])

        logging.info(f"P2Rank-determined grid center: X={center_x}, Y={center_y}, Z={center_z}")
        print(f"Using P2Rank-determined grid center: X={center_x} Å, Y={center_y} Å, Z={center_z} Å.")

        center_x += offset_x
        center_y += offset_y
        center_z += offset_z

        if offset_x != 0.0 or offset_y != 0.0 or offset_z != 0.0:
            logging.info(f"Applied offsets: ΔX={offset_x}, ΔY={offset_y}, ΔZ={offset_z}")
            print(f"Applied offsets: ΔX={offset_x} Å, ΔY={offset_y} Å, ΔZ={offset_z} Å.")
            logging.info(f"Shifted grid center: X={center_x}, Y={center_y}, Z={center_z}")
            print(f"Shifted grid center: X={center_x} Å, Y={center_y} Å, Z={center_z} Å.")

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
             exhaustiveness, energy_range, num_modes, seed, flex_pdbqt=None):
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
            '--num_modes', str(num_modes),
            '--seed', str(seed),
        ]

        if flex_pdbqt is not None and flex_pdbqt.strip():
            vina_command.extend(['--flex', flex_pdbqt])
            logging.info(f"Using flexible receptor: {flex_pdbqt}")

        logging.info(f"Vina command: {' '.join(vina_command)}")

        subprocess.run(vina_command, check=True)
        logging.info(f"Docking completed for ligand: {ligand_pdbqt}")

        affinities = []
        with open(output_pdbqt, 'r') as f:
            for line in f:
                if line.startswith('REMARK VINA RESULT:'):
                    parts = line.strip().split()
                    if len(parts) >= 6:
                        affinity = float(parts[3])
                        rmsd_lb = float(parts[4])
                        rmsd_ub = float(parts[5])
                        affinities.append((affinity, rmsd_lb, rmsd_ub))
        if affinities:
            table_output = 'mode |   affinity | rmsd l.b.| rmsd u.b.\n'
            table_output += '-----+------------+----------+----------\n'
            for idx, (affinity, rmsd_lb, rmsd_ub) in enumerate(affinities, start=1):
                table_output += f'{idx:>4}    {affinity:>10.4f}   {rmsd_lb:>8.3f}   {rmsd_ub:>8.3f}\n'
        else:
            table_output = "No docking results found."

        return table_output, affinities

    except subprocess.CalledProcessError as e:
        logging.error(f"Error in docking with Vina: {e}")
        print(f"Error in docking with Vina: {e}")
        raise

def add_axes(center=(0, 0, 0), length=10.0, radius=0.1, offset=(15, 15, 15)):
    x, y, z = center
    ox, oy, oz = offset
    new_center = (x + ox, y + oy, z + oz)

    axes = [
        CYLINDER, new_center[0], new_center[1], new_center[2], new_center[0] + length, new_center[1], new_center[2], radius,
        1.0, 0.0, 0.0,
        1.0, 0.0, 0.0,
        CYLINDER, new_center[0], new_center[1], new_center[2], new_center[0], new_center[1] + length, new_center[2], radius,
        0.0, 1.0, 0.0,
        0.0, 1.0, 0.0,
        CYLINDER, new_center[0], new_center[1], new_center[2], new_center[0], new_center[1], new_center[2] + length, radius,
        0.0, 0.0, 1.0,
        0.0, 0.0, 1.0
    ]

    cmd.load_cgo(axes, "axes")
    cmd.pseudoatom(object='axis_x_end', pos=(new_center[0] + length, new_center[1], new_center[2]))
    cmd.pseudoatom(object='axis_y_end', pos=(new_center[0], new_center[1] + length, new_center[2]))
    cmd.pseudoatom(object='axis_z_end', pos=(new_center[0], new_center[1], new_center[2] + length))
    cmd.label('axis_x_end', '"X"')
    cmd.label('axis_y_end', '"Y"')
    cmd.label('axis_z_end', '"Z"')

@logger_decorator
def generate_visualizations(receptor_pdbqt, output_pdbqt, output_folder, receptor_name, ligand_name, center, size, offset=(0, 0, 0)):
    try:
        cmd.load(receptor_pdbqt, 'receptor')
        cmd.load(output_pdbqt, 'ligand')
        cmd.h_add('receptor')
        cmd.h_add('ligand')

        cmd.color('cyan', 'receptor')
        cmd.color('red', 'ligand')
        cmd.show('cartoon', 'receptor')
        cmd.show('sticks', 'ligand')
        cmd.hide('sticks', '(elem H)')

        add_axes(center=center, length=10.0, radius=0.3, offset=offset)
        draw_docking_box(center, size)
        set_visualization_and_focus()

        image_path = os.path.join(output_folder, f"{receptor_name}_{ligand_name}_docking.png")
        cmd.ray(1920, 1080)
        cmd.png(image_path, width=1920, height=1080, dpi=300)
        logging.info(f"Visualization saved: {image_path}")

        session_path = os.path.join(output_folder, f"{receptor_name}_{ligand_name}_docking.pse")
        cmd.save(session_path)
        logging.info(f"PyMOL session saved: {session_path}")

        cmd.delete('all')
        
    except Exception as e:
        logging.error(f"Error in generating visualization: {e}")
        print(f"Error in generating visualization: {e}")
        raise

def draw_docking_box(center, size):
    x, y, z = center
    sx, sy, sz = size
    min_x = x - sx / 2
    max_x = x + sx / 2
    min_y = y - sy / 2
    max_y = y + sy / 2
    min_z = z - sz / 2
    max_z = z + sz / 2

    box = [
        LINEWIDTH, 1.0,
        COLOR, 1.0, 0.0, 1.0,
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
    cmd.bg_color('black')
    cmd.center('ligand')
    cmd.zoom('ligand', buffer=50)
    cmd.move('z', 120)
    cmd.turn('x', -5)
    cmd.turn('y', -5)
    cmd.color('skyblue', 'receptor')
    cmd.color('orange', 'ligand')
    cmd.color('white', 'axes')
    cmd.color('red', 'axis_x_end')
    cmd.color('green', 'axis_y_end')
    cmd.color('blue', 'axis_z_end')
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
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 150)
        options = drawer.drawOptions()
        options.padding = 0.1
        options.fixedFontSize = 11
        options.useFixedFontSize = True
        options.minFontSize = 6
        options.bondLineWidth = 2
        drawer.SetDrawOptions(options)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        svg = svg.replace('<?xml version=\'1.0\' encoding=\'utf-8\'?>\n', '')
        svg = svg.replace('xmlns:svg=', 'xmlns=')
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
        p2rank_csv = predictions_csv
        df_p2rank = pd.read_csv(p2rank_csv)
        df_p2rank.columns = df_p2rank.columns.str.strip()
        df_p2rank = df_p2rank.map(lambda x: x.strip() if isinstance(x, str) else x)
        df_p2rank['score'] = df_p2rank['score'].astype(float).map("{:.2f}".format)
        df_p2rank['probability'] = df_p2rank['probability'].astype(float).map("{:.2f}".format)

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
            hf.write('.probability { background-color: #ffecd9; }\n')
            hf.write('.docking-energy { background-color: #dfffe0; }\n')
            hf.write('td:nth-child(6), th:nth-child(6) {\n')
            hf.write('  max-width: 200px;\n')
            hf.write('  word-wrap: break-word;\n')
            hf.write('  white-space: normal;\n')
            hf.write('  text-align: center;\n')
            hf.write('}\n')
            hf.write('.p2rank-table td:nth-child(3), .p2rank-table th:nth-child(3) {\n')
            hf.write('  white-space: nowrap;\n')
            hf.write('  padding-left: 15px;\n')
            hf.write('  padding-right: 15px;\n')
            hf.write('}\n')
            hf.write('.p2rank-table td:nth-child(5), .p2rank-table th:nth-child(5) {\n')
            hf.write('  max-width: 400px;\n')
            hf.write('  word-wrap: break-word;\n')
            hf.write('  white-space: normal;\n')
            hf.write('  padding: 15px;\n')
            hf.write('  text-align: left;\n')
            hf.write('}\n')
            hf.write('</style>\n')
            hf.write('</head>\n')
            hf.write('<body>\n')

            receptor_pdbqt = os.path.basename(receptor_pdbqt)
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

            hf.write('<table>\n')
            hf.write(
                '<tr>'
                '<th>Number</th>'
                '<th>Compound Name</th>'
                '<th>Structure</th>'
                '<th>Docking Image</th>'
                '<th class="docking-energy">Affinity<br/>(kcal/mol)</th>'
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
                output_pdbqt_path = os.path.relpath(result['output_pdbqt'], os.path.dirname(html_file))
                link_text = os.path.basename(result['output_pdbqt'])
                pymol_session_path = os.path.relpath(result['pymol_session'], os.path.dirname(html_file))
                pymol_session_link = os.path.basename(result['pymol_session'])

                hf.write('<tr>\n')
                hf.write(f'<td>{idx}</td>\n')
                hf.write(f'<td>{name}</td>\n')
                hf.write(f'<td><img src="{image_path}" alt="{name}" width="400"/></td>\n')
                hf.write(
                    f'<td><a href="{docking_image_path}" target="_blank">'
                    f'<img src="{docking_image_path}" alt="Docking Image" style="max-width: 150px; max-height: 150px;"></a></td>\n'
                )
                hf.write(f'<td class="docking-energy">{affinity_str}</td>\n')
                hf.write(
                    f'<td><a href="{output_pdbqt_path}" download="{link_text}" type="application/octet-stream">'
                    f'{link_text}</a></td>\n'
                )
                hf.write(
                    f'<td><a href="{pymol_session_path}" download="{pymol_session_link}" '
                    f'type="application/octet-stream">{pymol_session_link}</a></td>\n'
                )
                hf.write('</tr>\n')
            hf.write('</table>\n')
            hf.write('</br>')

            results_file = f"{receptor_name}_results.txt"
            hf.write('<p style="text-align: center; margin-top: 20px;">\n')
            hf.write(
                f'<a href="{results_file}" target="_blank">'
                f'Detailed results for each compound (All docking poses). CLICK</a>\n'
            )
            hf.write('</p>\n')
            hf.write('</br></br>')

            p2rank_header = f'P2RANK: identified docking pockets for receptor with PDB code: {receptor_name}'
            hf.write(f'<h3 style="text-align: center; margin-top: 20px;">{p2rank_header}</h3>\n')

            hf.write('<table class="p2rank-table">\n')
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
            residues_csv_file = os.path.join('01_p2rank_output', f'{receptor_name}_fixed.pdb_residues.csv')
            hf.write('<p style="text-align: center; margin-top: 20px;">\n')
            hf.write(
                f'<a href="{residues_csv_file}" target="_blank">'
                f'Detailed information about individual amino acids </br>and their involvement in docking pockets. CLICK</a>\n'
            )
            hf.write('</p>\n')

            hf.write('</br></br>')

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

if __name__ == "__main__":
    main()
