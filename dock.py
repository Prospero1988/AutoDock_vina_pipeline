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
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Ścieżki do narzędzi
P2RANK_PATH = "/usr/local/bin/prank"
VINA_PATH = "/usr/local/bin/vina_1.2.5_linux_x86_64"
OBABEL_PATH = "/usr/bin/obabel"

# Sprawdzenie, czy narzędzia są dostępne
if not os.path.exists(P2RANK_PATH):
    raise FileNotFoundError(f"P2Rank executable not found at {P2RANK_PATH}")
if not os.path.exists(VINA_PATH):
    raise FileNotFoundError(f"AutoDock Vina executable not found at {VINA_PATH}")
if shutil.which(OBABEL_PATH) is None:
    raise EnvironmentError(f"Open Babel not found at {OBABEL_PATH}")

# Inicjalizacja PyMOL w trybie cichym bez GUI
pymol.finish_launching(['pymol', '-qc'])

folder_name = 'retinol_binding'
receptor_name = 'retinolP'
PDB_ID = '1RBP'
ligand_name = 'retinol'
pubchem_ID = 445354

# Tworzenie folderu o podanej nazwie
if not os.path.exists(folder_name):
    os.makedirs(folder_name, exist_ok=True)
    print(f"Folder '{folder_name}' został utworzony pomyślnie!")
else:
    print(f"Folder '{folder_name}' już istnieje.")

# Wyświetlanie wprowadzonych nazw
print(f"Nazwa receptora: {receptor_name}")
print(f"Nazwa liganda: {ligand_name}")

# Funkcja pobierająca SMILES z PubChem CID
def fetch_smiles_from_cid(cid):
    compound = pcp.Compound.from_cid(cid)
    return compound.isomeric_smiles

def generate_minimized_pdb(smiles, pdb_filename):
    # Generowanie molekuły z SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Nie udało się wygenerować molekuły z SMILES.")
    else:
        print("Molekuła została poprawnie wygenerowana z SMILES.")

    # Dodawanie atomów wodoru
    mol = Chem.AddHs(mol)

    # Generowanie współrzędnych 3D
    AllChem.EmbedMolecule(mol)

    # Próbuj najpierw z MMFF
    if AllChem.MMFFHasAllMoleculeParams(mol):
        print("Rozpoczęto minimalizację energii przy użyciu MMFF...")
        result = AllChem.MMFFOptimizeMolecule(mol)
        if result == 0:
            raise RuntimeError("Minimalizacja energii przy użyciu MMFF nie powiodła się.")
        print("Minimalizacja energii przy użyciu MMFF zakończona sukcesem.")
    else:
        print("MMFF nie jest obsługiwane dla tej molekuły. Próbuję z UFF...")
        # Jeśli MMFF nie działa, użyj UFF jako zapasowego
        result = AllChem.UFFOptimizeMolecule(mol)
        if result != 0:
            raise RuntimeError("Minimalizacja energii przy użyciu UFF nie powiodła się.")
        print("Minimalizacja energii przy użyciu UFF zakończona sukcesem.")

    # Sanitacja molekuły
    Chem.SanitizeMol(mol)

    # Zapis jako PDB
    Chem.MolToPDBFile(mol, pdb_filename)
    print(f"Zminimalizowana molekuła została zapisana jako {pdb_filename}")

cid = pubchem_ID
pdb_filename = f'{folder_name}/{ligand_name}.pdb'
smiles = fetch_smiles_from_cid(cid)
generate_minimized_pdb(smiles, pdb_filename)

def download_pdb(pdb_id, download_dir):
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)

    pdbl = PDBList()
    pdb_file_path = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb', pdir=download_dir)
    if not os.path.exists(pdb_file_path):
        raise FileNotFoundError(f"Failed to download PDB {pdb_id}.")

    return pdb_file_path

downloaded_pdb_path = download_pdb(PDB_ID, folder_name)

source = downloaded_pdb_path
destination = f'{folder_name}/{receptor_name}_dirty.pdb'
shutil.move(source, destination)

def fix_pdb(input_pdb, output_pdb, ph=7.0, chain_id='A'):
    fixer = PDBFixer(filename=input_pdb)
    
    # Usuń wszystkie łańcuchy poza wskazanym
    fixer.removeChains([chain.id for chain in fixer.topology.chains() if chain.id != chain_id])
    fixer.findMissingResidues()

    # Pobierz listę reszt z wybranego łańcucha
    chain_residues = list(fixer.topology.chains(chain_id).residues)
    print(f"Reszty w łańcuchu {chain_id}: {[res for res in fixer.topology.chains(chain_id).residues]}")

    # Usuń brakujące reszty na końcach łańcucha
    fixer.missingResidues = {
        key: val
        for key, val in fixer.missingResidues.items()
        if not (key[1] == 0 or key[1] == len(chain_residues) - 1)
    }

    # Znajdź i dodaj brakujące atomy oraz wodory
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)

    # Zapisz wynik do pliku
    with open(output_pdb, 'w') as outfile:
        PDBFile.writeFile(fixer.topology, fixer.positions, outfile)

    print(f"Naprawiony PDB zapisany jako {output_pdb}")



fixed_pdb = f'{folder_name}/{receptor_name}_fixed.pdb'
fix_pdb(f'{folder_name}/{receptor_name}_dirty.pdb', fixed_pdb, ph=7.4, chain_id='A')

def prepare_receptor(input_pdb, output_pdbqt):
    subprocess.run([OBABEL_PATH, "-i", "pdb", input_pdb, "-o", "pdbqt", "-O", output_pdbqt, "-xr"], check=True)
    print(f"Receptor przygotowany: {output_pdbqt}")

receptor_pdb = fixed_pdb
receptor_pdbqt = f"{folder_name}/{receptor_name}.pdbqt"
prepare_receptor(receptor_pdb, receptor_pdbqt)

ligand_pdb = f"{folder_name}/{ligand_name}.pdb"
ligand_pdbqt = f"{folder_name}/{ligand_name}.pdbqt"

def prepare_ligand(input_pdb, output_pdbqt):
    subprocess.run([OBABEL_PATH, "-i", "pdb", input_pdb, "-o", "pdbqt", "-O", output_pdbqt, "-h"], check=True)
    print(f"Ligand przygotowany: {output_pdbqt}")

prepare_ligand(ligand_pdb, ligand_pdbqt)

subprocess.run([P2RANK_PATH, 'predict', '-f', receptor_pdb])

receptor_base_name = os.path.basename(receptor_pdb).split('.')[0]
df = pd.read_csv(f'p2rank_2.4.2/test_output/predict_{receptor_base_name}/{receptor_base_name}.pdb_predictions.csv')
center_x = float(df['   center_x'].iloc[0])
center_y = float(df['   center_y'].iloc[0])
center_z = float(df['   center_z'].iloc[0])

pred = pd.read_csv(f'p2rank_2.4.2/test_output/predict_{receptor_base_name}/{receptor_base_name}.pdb_residues.csv')
pocket1 = pred[pred[' pocket'] == 1]
resi = '+'.join([str(i) for i in pocket1[' residue_label']])

cmd.load(receptor_pdb)
cmd.select('pocket1', f'resi {resi}')
cmd.show('cartoon')

alpha_carbons = []
model = cmd.get_model('pocket1 and name CA')
for atom in model.atom:
    alpha_carbons.append([atom.coord[0], atom.coord[1], atom.coord[2]])

alpha_carbons = np.array(alpha_carbons)
min_coords = np.min(alpha_carbons, axis=0)
max_coords = np.max(alpha_carbons, axis=0)
cube_size = max_coords - min_coords + 10
Size_x, Size_y, Size_z = cube_size

print(f"Center coordinates: x={center_x}, y={center_y}, z={center_z}")
print(f"Grid size: x={Size_x}, y={Size_y}, z={Size_z}")

output = f"{receptor_name}_{ligand_name}.pdbqt"
vina_command = [
    VINA_PATH,
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
subprocess.run(vina_command, check=True)

cmd.reinitialize()

cmd.load(f'{folder_name}/{output}')
cmd.load(receptor_pdb)

cmd.show('cartoon')
cmd.color('cyan', 'all')
cmd.color('red', f"{output.split('.')[0]}")
cmd.set('ray_trace_frames', 1)

output_image_path = f'{folder_name}/{receptor_name}_{ligand_name}_image.png'
cmd.png(output_image_path)

combined_pdb_path = f'{folder_name}/{receptor_name}_{ligand_name}_best.pdb'
cmd.save(combined_pdb_path, 'all', -1)

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
