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
- --exhaust: Specifies how thorough the search should be for the best binding poses.
                Higher values increase precision but require more computation time. (default: 20)
- --energy_range: Determines the range of energy scores (in kcal/mol) for poses to be considered. (default: 2)
- --keep_water: If specified, retains water in the receptor structure. By default, water is removed.

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
    parser = argparse.ArgumentParser(description="Automated docking script for multiple ligands.")
    parser.add_argument('--pdb_id', required=True, help='PDB ID of the receptor protein.')
    parser.add_argument('--ligands', required=True, help='Name of the SDF file containing ligands, located in ./ligands.')
    parser.add_argument('--tol', type=int, default=0, help='Tolerance in Ångströms to expand the docking pocket dimensions beyond those defined by P2Rank (default: 0).')
    parser.add_argument('--pckt', type=int, default=1, help='Pocket number to use from P2Rank predictions (default: 1).')
    parser.add_argument('--exhaust', type=int, default=20, help='Specifies how thorough the search should be for the best binding poses. Higher values increase precision but require more computation time (default: 20).')
    parser.add_argument('--energy_range', type=int, default=2, help='Determines the range of energy scores (in kcal/mol) for poses to be considered (default: 2).')
    parser.add_argument('--keep_water', action='store_true', help='If specified, retains water molecules in the receptor structure. By default, water is removed.')
    parser.add_argument('--keep_heterogens', action='store_true', help='If specified, retains heteroatoms (e.g., ions, cofactors) in the receptor PDB structure. By default, heteroatoms are removed.')

    args = parser.parse_args()
    pocket_number = args.pckt
    PDB_ID = args.pdb_id
    receptor_name = PDB_ID
    ligands_file = args.ligands
    tol = args.tol
    exhaustiveness = args.exhaust
    energy_range = args.energy_range
    keepWater = args.keep_water
    heterogens = args.keep_heterogens

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

    # Initialize list to collect ligand results
    ligand_results = []

    # Start processing
    try:
        # Download and prepare receptor
        downloaded_pdb_path = download_pdb(PDB_ID, receptor_folder)
        dirty_pdb = os.path.join(receptor_folder, f'{receptor_name}_dirty.pdb')
        shutil.move(downloaded_pdb_path, dirty_pdb)

        fixed_pdb = os.path.join(receptor_folder, f'{receptor_name}_fixed.pdb')
        fix_pdb(dirty_pdb, fixed_pdb, heterogens, keepWater, ph=7.4)

        receptor_pdbqt = os.path.join(receptor_folder, f"{receptor_name}.pdbqt")
        prepare_receptor(fixed_pdb, receptor_pdbqt)

        # Run P2Rank
        output_dir = os.path.join(receptor_folder, '01_p2rank_output')
        run_p2rank(fixed_pdb, output_dir)

        # Get docking box parameters
        center_x, center_y, center_z, Size_x, Size_y, Size_z, predictions_csv = get_docking_box(output_dir, fixed_pdb, tol, pocket_number)

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

                # Get the first affinity value
                if affinities:
                    affinity = affinities[0][0]
                else:
                    affinity = None

                # Collect ligand results
                ligand_results.append({'name': ligand_name, 'image': image_filename, 'affinity': affinity})

        # Generate HTML results file
        html_file = os.path.join(receptor_folder, f"{receptor_name}_results.html")
        generate_html_results(html_file, receptor_name, ligands_file, ligand_results, predictions_csv)

        print(f'Results HTML file generated: {html_file}')

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
def fix_pdb(input_pdb, output_pdb, heterogens, keepWater, ph=7.4):
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
        if heterogens:
            logging.info("Heterogens are retained in the receptor structure.")
            print("Heterogens are retained in the receptor structure.")
            # Do not call removeHeterogens() since we want to keep heterogens
            if not keepWater:
                # Remove only water molecules
                fixer.removeHeterogens(keepWater=False)
                logging.info("Water molecules removed from the receptor structure.")
                print("Water molecules are removed from the receptor structure.")
        else:
            # Remove heterogens; control water retention with keepWater
            fixer.removeHeterogens(keepWater=keepWater)
            logging.info(f"Heterogens removed; water molecules {'retained' if keepWater else 'removed'}.")
            print(f"Heterogens removed; water molecules {'retained' if keepWater else 'removed'} in the receptor structure.")

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

        return table_output, affinities  # Zmieniono, aby zwracać również listę affinities

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

def draw_molecule_to_file(mol, image_filename, width=400):
    try:
        # Oblicz bounding box molekuły
        conf = mol.GetConformer()
        coords = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
        x_coords = [coord.x for coord in coords]
        y_coords = [coord.y for coord in coords]

        # Rozpiętość w osiach X i Y
        x_range = max(x_coords) - min(x_coords)
        y_range = max(y_coords) - min(y_coords)

        # Obliczenie wysokości proporcjonalnie do szerokości
        if x_range > 0:
            aspect_ratio = y_range / x_range
        else:
            aspect_ratio = 1.0  # Domyślna proporcja, gdyby coś poszło nie tak
        height = int(width * aspect_ratio)

        # Przygotowanie rysownika SVG
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        options = drawer.drawOptions()
        options.padding = 0.1  # Mały margines
        options.fixedFontSize = 11  # Stały rozmiar czcionki
        options.useFixedFontSize = True  # Wymuszenie stałego rozmiaru czcionki
        options.minFontSize = 6  # Minimalny rozmiar czcionki
        options.bondLineWidth = 2  # Grubsze linie wiązań
        drawer.SetDrawOptions(options)

        # Rysowanie molekuły
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()

        # Usunięcie deklaracji XML i poprawa przestrzeni nazw
        svg = svg.replace('<?xml version=\'1.0\' encoding=\'utf-8\'?>\n', '')
        svg = svg.replace('xmlns:svg=', 'xmlns=')

        # Zapisanie pliku SVG
        with open(image_filename, 'w') as f:
            f.write(svg)
        logging.info(f"Ligand image saved: {image_filename}")
    except Exception as e:
        logging.error(f"Error in drawing molecule image: {e}")
        print(f"Error in drawing molecule image: {e}")
        raise

def generate_html_results(html_file, receptor_name, ligands_file, ligand_results, predictions_csv):
    try:
        # Wczytaj dane z P2RANK i wyczyść białe znaki
        p2rank_csv = predictions_csv
        df_p2rank = pd.read_csv(p2rank_csv)
        df_p2rank.columns = df_p2rank.columns.str.strip()
        df_p2rank = df_p2rank.map(lambda x: x.strip() if isinstance(x, str) else x)
        df_p2rank['score'] = df_p2rank['score'].astype(float).map("{:.2f}".format)

        with open(html_file, 'w', encoding='utf-8') as hf:
            hf.write('<html>\n')
            hf.write('<head>\n')
            hf.write('<title>Wyniki dokowania</title>\n')
            hf.write('<style>\n')
            hf.write('body { background-color: white; font-family: Arial, sans-serif; }\n')
            hf.write('table { border-collapse: collapse; margin: auto; }\n')
            hf.write('th, td { border: 1px solid black; padding: 5px; text-align: center; vertical-align: middle; }\n')
            hf.write('th { background-color: #f2f2f2; }\n')
            hf.write('td:nth-child(3), th:nth-child(3) { width: 400px; } /* Kolumna z obrazkami */\n')
            hf.write('img { display: block; margin: auto; }\n')
            hf.write('.probability { background-color: #ffecd9; } /* Pastelowy pomarańczowy */\n')
            # Dodano style dla drugiej tabelki, aby dostosować szerokość kolumny "score"
            
            hf.write('.p2rank-table td:nth-child(3), .p2rank-table th:nth-child(3) { max-width: 100px; white-space: nowrap; text-align: center; }\n')
            # Styl dla drugiej tabelki (dostosowanie kolumny residue_ids)
            hf.write('.p2rank-table td:nth-child(5) {\n')
            hf.write('  max-width: 300px;\n')
            hf.write('  word-wrap: break-word;\n')  # Zawijanie tekstu
            hf.write('  text-align: left;\n')  # Wyrównanie do lewej
            hf.write('  vertical-align: top;\n')  # Wyrównanie do góry
            hf.write('}\n')
           
            hf.write('</style>\n')
            hf.write('</head>\n')
            hf.write('<body>\n')
                        
            # Nagłówek dla tabelki z wynikami dokowania
            header_text = f'Wyniki dokowania do receptora o kodzie PDB: {receptor_name} </br>przy pomocy struktur zapisanych w pliku: {ligands_file}</br>'
            hf.write(f'<h2 style="text-align: center;">{header_text}</h2>\n')

            # Pierwsza tabela: Wyniki dokowania
            hf.write('<table>\n')
            hf.write('<tr><th>Numer</th><th>Nazwa związku</th><th>Struktura</th><th>Energia dokowania<br/>(kcal/mol)</th></tr>\n')
            for idx, result in enumerate(ligand_results, start=1):
                name = result['name']
                image_path = os.path.relpath(result['image'], os.path.dirname(html_file))
                affinity = result['affinity']
                if affinity is not None:
                    affinity_str = f"{affinity:.2f}"
                else:
                    affinity_str = 'N/A'
                hf.write('<tr>\n')
                hf.write(f'<td>{idx}</td>\n')
                hf.write(f'<td>{name}</td>\n')
                hf.write(f'<td><img src="{image_path}" alt="{name}" width="400"/></td>\n')
                hf.write(f'<td>{affinity_str}</td>\n')
                hf.write('</tr>\n')
            hf.write('</table>\n')
            hf.write('</br>')
            
            # Link do szczegółowych rezultatów
            results_file = f"{receptor_name}_results.txt"
            hf.write('<p style="text-align: center; margin-top: 20px;">\n')
            hf.write(f'<a href="{results_file}" target="_blank">Szczegółowe rezultaty dla każdego związku (Wszystkie pozy dokowania). KLIKNIJ</a>\n')
            hf.write('</p>\n')
            hf.write('</br></br>')
            
            # Nagłówek dla tabelki z danymi P2RANK
            p2rank_header = f'P2RANK: wyznaczone kieszenie dokowania dla receptora o kodzie PDB: {receptor_name}'
            hf.write(f'<h3 style="text-align: center; margin-top: 20px;">{p2rank_header}</h3>\n')

            # Druga tabela: Dane P2RANK
            hf.write('<table class="p2rank-table">\n')  # Dodano klasę "p2rank-table"
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
            # Link do szczegółowych rezultatów
            residues_csv_file = os.path.join('01_p2rank_output', f'{receptor_name}_fixed.pdb_residues.csv')
            hf.write('<p style="text-align: center; margin-top: 20px;">\n')
            hf.write(f'<a href="{residues_csv_file}" target="_blank">Szczegółowe informacje o poszczególnych aminokwasach i ich udziale w kieszeniach dokowania. KLIKNIJ</a>\n')
            hf.write('</p>\n')
            
            hf.write('</br></br>')
            
            # Informacje o autorze na końcu
            hf.write('<p style="font-size: small; text-align: center; margin-top: 20px;">\n')
            hf.write('System dokowania oparty o <b>AutoDock Vina v.1.2.5</b> oraz <b>P2RANK v.2.4.2</b><br/>\n')
            hf.write('<b>Autor:</b> Arkadiusz Leniak <b>mail:</b> arkadiusz.leniak@gmail.com<br/>\n')
            hf.write('<b>github:</b> <a href="https://github.com/Prospero1988">https://github.com/Prospero1988</a>\n')
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
