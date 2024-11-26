from Bio import PDB
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem
import os
from Bio.PDB import PDBList
import pandas as pd
import numpy as np
import pymol
from pymol import cmd
import shutil
import subprocess

# Initialize PyMOL in quiet mode without GUI
pymol.finish_launching(['pymol', '-qc'])

folder_name = 'retinol_binding'
receptor_name = 'retinolP'
PDB_ID = '1RBP'
ligand_name = 'retinol'
pubchem_ID = 445354

# Create the folder with the specified name
if not os.path.exists(folder_name):
    os.mkdir(folder_name)
    print(f"Folder '{folder_name}' created successfully!")
else:
    print(f"Folder '{folder_name}' already exists.")

# Displaying the entered names
print(f"Receptor Name: {receptor_name}")
print(f"Ligand Name: {ligand_name}")

# Function to fetch SMILES from PubChem CID
def fetch_smiles_from_cid(cid):
    compound = pcp.Compound.from_cid(cid)
    return compound.isomeric_smiles

def generate_minimized_pdb(smiles, pdb_filename):
    # Generate molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string provided.")

    # Add hydrogen atoms
    mol = Chem.AddHs(mol)

    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol)

    # Minimize energy
    AllChem.UFFOptimizeMolecule(mol)

    # Sanitize molecule
    Chem.SanitizeMol(mol)

    # Save as PDB
    Chem.MolToPDBFile(mol, pdb_filename)

    print(f"Minimized molecule saved as {pdb_filename}")

# Generate minimized ligand PDB
cid = pubchem_ID
pdb_filename = f'{folder_name}/{ligand_name}.pdb'
smiles = fetch_smiles_from_cid(cid)
generate_minimized_pdb(smiles, pdb_filename)

# Download and Pre-process Receptor
def download_pdb(pdb_id, download_dir):
    # Create the download directory if it doesn't exist
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)

    # Initialize PDB downloader
    pdbl = PDBList()

    # Download the PDB file
    pdb_file_path = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb', pdir=download_dir)

    return pdb_file_path

downloaded_pdb_path = download_pdb(PDB_ID, folder_name)

# Move and rename the downloaded PDB file
source = downloaded_pdb_path
destination = f'{folder_name}/{receptor_name}_dirty.pdb'
shutil.move(source, destination)

# Remove HETATM from PDB file
def remove_hetatm(input_pdb, output_pdb):
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if not line.startswith('HETATM'):
                outfile.write(line)

remove_hetatm(f'{folder_name}/{receptor_name}_dirty.pdb', f'{folder_name}/{receptor_name}_cleaner.pdb')

# Keep only chain A from PDB file
def filter_chain_a(input_pdb, output_pdb):
    # Initialize PDB parser and writer
    parser = PDB.PDBParser(QUIET=True)
    writer = PDB.PDBIO()

    # Parse the structure from the input PDB file
    structure = parser.get_structure('structure', input_pdb)

    # Create a new structure only containing chain A
    chain_a_model = PDB.Structure.Structure('filtered')
    for model in structure:
        new_model = PDB.Model.Model(model.id)
        chain_a = [chain for chain in model if chain.id == 'A']
        if chain_a:
            new_model.add(chain_a[0])
        chain_a_model.add(new_model)

    # Save the filtered structure to the output file
    writer.set_structure(chain_a_model)
    writer.save(output_pdb)
    print(f"Filtered PDB saved to {output_pdb}")

# Filter and save
filter_chain_a(f'{folder_name}/{receptor_name}_cleaner.pdb', f'{folder_name}/{receptor_name}.pdb')

# Define Box (p2rank)
if not os.path.exists('p2rank_2.4.2'):
    subprocess.run(['wget', 'https://github.com/rdk/p2rank/releases/download/2.4.2/p2rank_2.4.2.tar.gz'])
    subprocess.run(['tar', '-xzf', 'p2rank_2.4.2.tar.gz'])

subprocess.run(['p2rank_2.4.2/prank', 'predict', '-f', f'{folder_name}/{receptor_name}.pdb'])

df = pd.read_csv(f'p2rank_2.4.2/test_output/predict_{receptor_name}/{receptor_name}.pdb_predictions.csv')
center_x = float(df['   center_x'].iloc[0])
center_y = float(df['   center_y'].iloc[0])
center_z = float(df['   center_z'].iloc[0])

pred = pd.read_csv(f'p2rank_2.4.2/test_output/predict_{receptor_name}/{receptor_name}.pdb_residues.csv')
pocket1 = pred[pred[' pocket'] == 1]
resi = '+'.join([str(i) for i in pocket1[' residue_label']])

# Load receptor structure into PyMOL
cmd.load(f'{folder_name}/{receptor_name}.pdb')
cmd.select('pocket1', f'resi {resi}')
cmd.show('cartoon')

# Get alpha carbon coordinates
alpha_carbons = []
model = cmd.get_model('pocket1 and name CA')
for atom in model.atom:
    alpha_carbons.append([atom.coord[0], atom.coord[1], atom.coord[2]])

# Convert to numpy array for easier calculations
alpha_carbons = np.array(alpha_carbons)

# Find the min and max coordinates along each axis
min_coords = np.min(alpha_carbons, axis=0)
max_coords = np.max(alpha_carbons, axis=0)

# Calculate the dimensions of the bounding cube
cube_size = max_coords - min_coords
Size_x, Size_y, Size_z = cube_size

print(center_x, center_y, center_z, Size_x, Size_y, Size_z)

# Docking
def convert_pdb_to_pdbqt_receptor(input_pdb, output_pdbqt):
    subprocess.run(['obabel', '-i', 'pdb', input_pdb, '-o', 'pdbqt', '-O', output_pdbqt, '-xr', '-xn', '-xp'])

def prepare_ligand(input_pdb, output_pdb, ph):
    subprocess.run(['obabel', '-i', 'pdb', input_pdb, '-o', 'pdb', '-O', output_pdb, '--ph', str(ph), '-h'])

def minimize_energy(input_pdb, output_pdb, forcefield="MMFF94", steps=2500):
    subprocess.run(['obabel', input_pdb, '-O', output_pdb, '--minimize', '--ff', forcefield, '--steps', str(steps)])

def convert_pdb_to_pdbqt_ligand(input_pdb, output_pdbqt):
    subprocess.run(['obabel', '-i', 'pdb', input_pdb, '-o', 'pdbqt', '-O', output_pdbqt, '-h'])

# Define pH value, file paths, and other parameters
pH = 7.4
receptor_pdb = f"{folder_name}/{receptor_name}.pdb"
ligand_pdb = f"{folder_name}/{ligand_name}.pdb"
receptor_pdbqt = f"{folder_name}/{receptor_name}.pdbqt"
ligand_pdbqt = f"{folder_name}/{ligand_name}.pdbqt"

# Convert receptor PDB to PDBQT
convert_pdb_to_pdbqt_receptor(receptor_pdb, receptor_pdbqt)

# Convert minimized ligand to PDBQT
convert_pdb_to_pdbqt_ligand(ligand_pdb, ligand_pdbqt)

# Run AutoDock Vina
output = f"{receptor_name}_{ligand_name}.pdbqt"
vina_command = [
    'vina_1.2.5_linux_x86_64',
    '--receptor', receptor_pdbqt,
    '--ligand', ligand_pdbqt,
    '--out', f'{folder_name}/{output}',
    '--center_x', str(center_x),
    '--center_y', str(center_y),
    '--center_z', str(center_z),
    '--size_x', str(Size_x),
    '--size_y', str(Size_y),
    '--size_z', str(Size_z)
]
subprocess.run(vina_command)

# PyMol Visualization
cmd.reinitialize()

# Load the files
cmd.load(f'{folder_name}/{output}')
cmd.load(receptor_pdb)

# Display the structure
cmd.show('cartoon')

# Additional settings (optional)
cmd.color('cyan', 'all')  # Color the structure cyan
cmd.color('red', f"{output.split('.')[0]}")
cmd.set('ray_trace_frames', 1)  # Enable ray tracing for better quality

# Save an image of the structure
output_image_path = f'{folder_name}/{receptor_name}_{ligand_name}_image.png'
cmd.png(output_image_path)

# Save the combined structure as a PDB file
combined_pdb_path = f'{folder_name}/{receptor_name}_{ligand_name}_best.pdb'
cmd.save(combined_pdb_path, 'all', -1)

# Alignment
cmd.reinitialize()
cmd.load(f'{folder_name}/{receptor_name}_{ligand_name}_best.pdb', 'docking')
cmd.color('cyan', 'docking')
cmd.load(f'{folder_name}/{receptor_name}_dirty.pdb', 'RealStructure')
cmd.color('green', 'RealStructure')
cmd.align('docking', 'RealStructure')
cmd.show('cartoon')
output_image_path = f'{folder_name}/alignment.png'
cmd.png(output_image_path)
combined_pdb_path = f'{folder_name}/aligned.pdb'
cmd.save(combined_pdb_path, 'all', -1)
