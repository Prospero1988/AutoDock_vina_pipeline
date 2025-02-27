#!/usr/bin/env python3
"""
MANUAL DOCKING SCRIPT - Extended
-------------------------------------
This script performs docking of multiple ligands onto multiple receptors,
where each receptor is already a ready-to-dock .pdbqt file.

Folders used:
- ./receptors/   : contains *.pdbqt files (one per receptor).
- ./ligands/     : contains exactly one .sdf or .mol2 file with all ligands.
- ./parameters/  : contains CSV files, each named exactly like the receptor (no extension),
                   but with .csv extension. For example:
                     Receptor: "1ABC.pdbqt"
                     Params:   "1ABC.csv"

Each CSV has parameters needed for docking, for example:

parameter,value
exhaust,40
energy_range,1
num_modes,5
seed,1988
grid_center_x,107
grid_center_y,114
grid_center_z,143
grid_size_x,17
grid_size_y,20
grid_size_z,35

Workflow:
  1. Finds the ligand file in ./ligands (SDF or MOL2).
  2. Reads all receptor .pdbqt files in ./receptors.
  3. For each receptor, reads the matching CSV in ./parameters.
  4. For each ligand in the ligand file, converts it to .pdbqt and runs docking with AutoDock Vina,
     using the manually-defined grid parameters from the CSV.
  5. Generates:
     - a CSV with summarized docking results ( name, affinity, smiles ),
     - a text file ( "receptor_name_results.txt" ) with all detailed docking poses,
     - an HTML report with large images, link to .pdbqt, etc. (no P2Rank references).

Requires:
- RDKit
- Open Babel
- PyMOL + PyMol libraries
- AutoDock Vina
"""

import os
import re
import csv
import shutil
import logging
import subprocess
import tempfile
from pathlib import Path

import pymol
from pymol import cmd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

from pymol.cgo import CYLINDER
from pymol.cgo import BEGIN, LINES, VERTEX, END, COLOR, LINEWIDTH

# ----------------------------------------------
# User-configurable paths to external tools:
# ----------------------------------------------
VINA_PATH = "/usr/local/bin/vina_1.2.5_linux_x86_64"
OBABEL_PATH = "/usr/bin/obabel"

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')
pymol.finish_launching(['pymol', '-qc'])  # Start PyMOL in quiet mode (no GUI)

# -------------------
# Decorator for logging
# -------------------
def logger_decorator(func):
    def wrapper(*args, **kwargs):
        logging.info(f"Running {func.__name__}")
        return func(*args, **kwargs)
    return wrapper

@logger_decorator
def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # 1) Identify the single ligand file in ./ligands (SDF or MOL2)
    ligands_dir = os.path.join(script_dir, "ligands")
    if not os.path.exists(ligands_dir):
        raise FileNotFoundError(f"'ligands' folder not found at: {ligands_dir}")

    ligand_file = find_ligand_file(ligands_dir)
    ligand_molecules = read_ligands(ligand_file)
    logging.info(f"Found {len(ligand_molecules)} ligands in: {ligand_file}")

    # 2) For each receptor in ./receptors
    receptors_dir = os.path.join(script_dir, "receptors")
    if not os.path.exists(receptors_dir):
        raise FileNotFoundError(f"'receptors' folder not found at: {receptors_dir}")

    processed_receptors = []

    # Collect all .pdbqt files
    receptor_files = [f for f in os.listdir(receptors_dir) if f.lower().endswith(".pdbqt")]
    if not receptor_files:
        raise FileNotFoundError(f"No .pdbqt receptor files found in {receptors_dir}")

    parameters_dir = os.path.join(script_dir, "parameters")
    if not os.path.exists(parameters_dir):
        raise FileNotFoundError(f"'parameters' folder not found at: {parameters_dir}")

    # Loop over each receptor
    for receptor_file in receptor_files:
        receptor_name = os.path.splitext(receptor_file)[0]  # e.g. "1ABC" if file is "1ABC.pdbqt"
        processed_receptors.append(receptor_name)
        receptor_path = os.path.join(receptors_dir, receptor_file)

        # 3) Matching CSV in ./parameters
        param_csv_path = os.path.join(parameters_dir, f"{receptor_name}.csv")
        if not os.path.exists(param_csv_path):
            logging.warning(f"No parameter CSV found for receptor '{receptor_name}'. Skipping.")
            continue

        # Read the parameters from CSV
        param_dict = read_parameters(param_csv_path)
        logging.info(f"Read parameter file: {param_csv_path}")

        # Create a folder for this receptor's results
        receptor_folder = os.path.join(script_dir, receptor_name)
        os.makedirs(receptor_folder, exist_ok=True)

        # Folder for ligand results
        ligands_results_folder = os.path.join(receptor_folder, "ligands_results")
        os.makedirs(ligands_results_folder, exist_ok=True)

        # We'll collect results in a list for final HTML
        docking_results = []

        # Also create a text file that accumulates docking poses
        results_text_file = os.path.join(receptor_folder, f"{receptor_name}_results.txt")
        with open(results_text_file, "w", encoding="utf-8") as rf:
            rf.write(f"# Docking results for {receptor_name}\n\n")

        logging.info(f"Starting docking for receptor: {receptor_name}")
        # 4) For each ligand, dock with Vina
        for idx, mol in enumerate(ligand_molecules, 1):
            if mol is None:
                continue

            # 4A) Assign or sanitize ligand name
            lig_name = get_ligand_name(mol, idx)
            mol.SetProp("_Name", lig_name)

            # Generate SMILES for CSV
            smiles_str = Chem.MolToSmiles(mol)

            # Make a subfolder for this ligand
            ligand_folder = os.path.join(ligands_results_folder, lig_name)
            os.makedirs(ligand_folder, exist_ok=True)

            # Write ligand to .pdb
            ligand_pdb = os.path.join(ligand_folder, f"{lig_name}.pdb")
            write_ligand_to_pdb(mol, ligand_pdb)

            # Generate 2D structure (SVG)
            image_filename = os.path.join(ligand_folder, f"{lig_name}.svg")
            draw_molecule_to_file(mol, image_filename)

            # Convert .pdb to .pdbqt for docking
            ligand_pdbqt = os.path.join(ligand_folder, f"{lig_name}.pdbqt")
            prepare_ligand_for_docking(ligand_pdb, ligand_pdbqt)

            # Output file from Vina
            output_pdbqt = os.path.join(ligand_folder, f"{receptor_name}_{lig_name}_docked.pdbqt")

            # 4B) Perform docking
            table_output, affinities = run_vina(
                receptor_path,
                ligand_pdbqt,
                output_pdbqt,
                param_dict
            )

            # Write all docking poses to the text file
            with open(results_text_file, "a", encoding="utf-8") as rf:
                rf.write(f"Ligand: {lig_name}\n")
                rf.write(table_output)
                rf.write("\n\n")

            # Extract the best affinity if found
            best_affinity = None
            if affinities:
                best_affinity = affinities[0][0]  # first is the best (lowest energy)

            # 4C) Generate PyMOL visualization
            docking_image = os.path.join(ligand_folder, f"{receptor_name}_{lig_name}.png")
            pymol_session = os.path.join(ligand_folder, f"{receptor_name}_{lig_name}.pse")

            try:
                generate_visualization(
                    receptor_path,
                    output_pdbqt,
                    docking_image,
                    pymol_session,
                    param_dict
                )
            except Exception as e:
                logging.warning(f"Visualization error for ligand {lig_name}: {e}")

            # Store info in docking_results for the HTML
            docking_results.append({
                "name": lig_name,
                "affinity": best_affinity,
                "smiles": smiles_str,
                "image": image_filename,       # 2D structure
                "docking_image": docking_image,
                "output_pdbqt": output_pdbqt,
                "pymol_session": pymol_session
            })

        # 5) Create CSV summary ( name, affinity, smiles )
        csv_results_file = os.path.join(receptor_folder, f"{receptor_name}_results_in_CSV.csv")
        try:
            # Sort by best_affinity ascending (None last)
            docking_sorted = sorted(docking_results, key=lambda x: (x['affinity'] is None, x['affinity']))

            with open(csv_results_file, 'w', newline='', encoding='utf-8') as csvfile:
                fieldnames = ['name', 'affinity', 'smiles']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                for item in docking_sorted:
                    aff_val = item['affinity']
                    aff_str = f"{aff_val:.2f}" if aff_val is not None else ""
                    row = {
                        'name': item['name'],
                        'affinity': aff_str,
                        'smiles': item['smiles'] if item['smiles'] else ""
                    }
                    writer.writerow(row)

            logging.info(f"CSV summary created: {csv_results_file}")
        except Exception as e:
            logging.error(f"Error writing summary CSV: {e}")

        # 6) Generate an HTML report
        html_path = os.path.join(receptor_folder, f"{receptor_name}_results.html")
        generate_html_results(
            html_path,
            receptor_name,
            ligand_file,     # original ligand file name
            docking_results,
            receptor_path    # link to receptor .pdbqt
        )
        logging.info(f"Docking finished for receptor: {receptor_name} - results in {html_path}")

        # 7) Copy all output_pdbqt files to "03_ligands_PDBQT" folder
        ligands_pdbqt_folder = os.path.join(receptor_folder, "03_ligands_PDBQT")
        os.makedirs(ligands_pdbqt_folder, exist_ok=True)

        copied_count = 0
        for item in docking_results:
            src_pdbqt = item["output_pdbqt"]
            if os.path.isfile(src_pdbqt):
                shutil.copy(src_pdbqt, ligands_pdbqt_folder)
                copied_count += 1

        logging.info(
            f"Copied {copied_count} .pdbqt files to: {ligands_pdbqt_folder}"
        )

    # After processing all receptors, create "receptors_manual.csv"
    csv_path = os.path.join(receptors_dir, "receptors_manual.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["receptors"])
        for name in processed_receptors:
            writer.writerow([name])

    print(f"Created '{csv_path}' with {len(processed_receptors)} receptors listed.")


@logger_decorator
def find_ligand_file(ligands_dir):
    """Find exactly one .sdf or .mol2 file in the ligands folder."""
    possible_files = []
    for fname in os.listdir(ligands_dir):
        lower = fname.lower()
        if lower.endswith(".sdf") or lower.endswith(".mol2"):
            possible_files.append(os.path.join(ligands_dir, fname))

    if not possible_files:
        raise FileNotFoundError("No .sdf or .mol2 files found in ./ligands folder.")
    if len(possible_files) > 1:
        logging.warning("Multiple ligand files found. Will use the first one: %s", possible_files[0])
    return possible_files[0]


@logger_decorator
def read_ligands(ligand_file):
    """Read molecules from SDF or MOL2 using RDKit and (if needed) Open Babel conversion."""
    if ligand_file.lower().endswith(".sdf"):
        suppl = Chem.SDMolSupplier(ligand_file)
        if not suppl:
            raise ValueError(f"Could not read molecules from {ligand_file}")
        return [m for m in suppl if m is not None]

    elif ligand_file.lower().endswith(".mol2"):
        with tempfile.NamedTemporaryFile(delete=False, suffix='.sdf') as tmp_sdf:
            tmp_sdf_path = tmp_sdf.name
        try:
            cmd_line = [
                OBABEL_PATH,
                "-i", "mol2", ligand_file,
                "-o", "sdf", "-O", tmp_sdf_path,
                "--aromatic"
            ]
            subprocess.run(cmd_line, check=True)
            suppl = Chem.SDMolSupplier(tmp_sdf_path)
            if not suppl:
                raise ValueError(f"Could not read molecules from converted SDF: {tmp_sdf_path}")
            return [m for m in suppl if m is not None]
        finally:
            if os.path.exists(tmp_sdf_path):
                os.remove(tmp_sdf_path)
    else:
        raise ValueError("Ligand file must be .sdf or .mol2.")


@logger_decorator
def get_ligand_name(mol, idx):
    """Get or create a ligand name from RDKit property _Name."""
    if mol.HasProp("_Name"):
        raw_name = mol.GetProp("_Name").strip()
        if raw_name:
            return sanitize_ligand_name(raw_name)
    return f"ligand_{idx:03d}"


@logger_decorator
def sanitize_ligand_name(name):
    """Replace invalid chars with underscores."""
    valid_pattern = re.compile(r'[^A-Za-z0-9_-]')
    sanitized = valid_pattern.sub('_', name)
    sanitized = sanitized.strip('_')
    if not sanitized:
        sanitized = "ligand_unnamed"
    return sanitized


@logger_decorator
def write_ligand_to_pdb(mol, pdb_path):
    """Save RDKit molecule to a PDB file (with 3D coords)."""
    try:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)
        Chem.MolToPDBFile(mol, pdb_path)
        logging.info(f"Saved ligand to PDB: {pdb_path}")
    except Exception as e:
        logging.error(f"Error writing ligand to PDB: {e}")
        raise


@logger_decorator
def draw_molecule_to_file(mol, image_filename):
    """Generate a 2D depiction (SVG) of the RDKit molecule."""
    try:
        AllChem.Compute2DCoords(mol)
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
        logging.info(f"2D structure saved as: {image_filename}")

    except Exception as e:
        logging.error(f"Error in generating 2D structure: {e}")


@logger_decorator
def prepare_ligand_for_docking(ligand_pdb, ligand_pdbqt):
    """Use Open Babel to convert ligand PDB -> PDBQT (with hydrogens)."""
    cmd_line = [
        OBABEL_PATH,
        "-i", "pdb", ligand_pdb,
        "-o", "pdbqt",
        "-O", ligand_pdbqt,
        "-h"
    ]
    try:
        subprocess.run(cmd_line, check=True)
        logging.info(f"Converted {ligand_pdb} -> {ligand_pdbqt}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Open Babel conversion error: {e}")
        raise


@logger_decorator
def read_parameters(param_csv_path):
    """Reads param CSV with columns: parameter,value."""
    param_dict = {}
    with open(param_csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            p = row["parameter"].strip()
            v = row["value"].strip()
            param_dict[p] = v
    return param_dict


@logger_decorator
def run_vina(receptor_pdbqt, ligand_pdbqt, output_pdbqt, params_dict):
    """
    AutoDock Vina docking using values from param_dict:
      exhaust, energy_range, num_modes, seed, grid_center_x, etc.
    """
    exhaustiveness = params_dict.get("exhaust", "8")
    energy_range  = params_dict.get("energy_range", "3")
    num_modes     = params_dict.get("num_modes", "9")
    seed          = params_dict.get("seed", "1988")

    center_x = float(params_dict.get("grid_center_x", "0"))
    center_y = float(params_dict.get("grid_center_y", "0"))
    center_z = float(params_dict.get("grid_center_z", "0"))

    size_x = float(params_dict.get("grid_size_x", "10"))
    size_y = float(params_dict.get("grid_size_y", "10"))
    size_z = float(params_dict.get("grid_size_z", "10"))

    vina_cmd = [
        VINA_PATH,
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--out", output_pdbqt,
        "--center_x", str(center_x),
        "--center_y", str(center_y),
        "--center_z", str(center_z),
        "--size_x", str(size_x),
        "--size_y", str(size_y),
        "--size_z", str(size_z),
        "--exhaustiveness", exhaustiveness,
        "--energy_range", energy_range,
        "--num_modes", num_modes,
        "--seed", str(seed)
    ]

    try:
        subprocess.run(vina_cmd, check=True)
        logging.info(f"Docking complete. Output: {output_pdbqt}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Vina docking error: {e}")
        raise

    # Parse Vina results
    table_output = ""
    affinities = []
    with open(output_pdbqt, 'r') as f:
        for line in f:
            if line.startswith('REMARK VINA RESULT:'):
                parts = line.strip().split()
                if len(parts) >= 6:
                    aff = float(parts[3])
                    rmsd_lb = float(parts[4])
                    rmsd_ub = float(parts[5])
                    affinities.append((aff, rmsd_lb, rmsd_ub))

    if affinities:
        table_output += "mode |   affinity | rmsd l.b.| rmsd u.b.\n"
        table_output += "-----+------------+----------+----------\n"
        for i, (aff, lb, ub) in enumerate(affinities, start=1):
            table_output += f"{i:>4}    {aff:>10.4f}   {lb:>8.3f}   {ub:>8.3f}\n"
    else:
        table_output = "No docking results found."

    return table_output, affinities


@logger_decorator
def add_axes(center=(0, 0, 0), length=10.0, radius=0.1):
    """
    Draw X/Y/Z axes in PyMOL, color-coded, with label pseudoatoms at the ends.
    """
    axes = [
        CYLINDER, center[0], center[1], center[2],
                   center[0] + length, center[1], center[2],
                   radius,
                   1.0, 0.0, 0.0,  # start color: red
                   1.0, 0.0, 0.0,  # end color: red

        CYLINDER, center[0], center[1], center[2],
                   center[0], center[1] + length, center[2],
                   radius,
                   0.0, 1.0, 0.0,  # start color: green
                   0.0, 1.0, 0.0,  # end color: green

        CYLINDER, center[0], center[1], center[2],
                   center[0], center[1], center[2] + length,
                   radius,
                   0.0, 0.0, 1.0,  # start color: blue
                   0.0, 0.0, 1.0   # end color: blue
    ]

    cmd.load_cgo(axes, "axes")
    # label pseudoatoms
    cmd.pseudoatom(object='axis_x_end', pos=(center[0] + length, center[1], center[2]))
    cmd.pseudoatom(object='axis_y_end', pos=(center[0], center[1] + length, center[2]))
    cmd.pseudoatom(object='axis_z_end', pos=(center[0], center[1], center[2] + length))
    cmd.label('axis_x_end', '"X"')
    cmd.label('axis_y_end', '"Y"')
    cmd.label('axis_z_end', '"Z"')


@logger_decorator
def generate_visualization(
    receptor_pdbqt,
    docked_ligand_pdbqt,
    png_output,
    pse_output,
    params_dict
):
    """
    PyMOL-based visualization with the X/Y/Z axes centered at the grid box center,
    plus the bounding box, plus saving .png and .pse.
    """
    cmd.reinitialize()
    cmd.load(receptor_pdbqt, "receptor")
    cmd.load(docked_ligand_pdbqt, "ligand")

    cmd.hide("everything", "all")
    cmd.show("cartoon", "receptor")
    cmd.show("sticks", "ligand")
    cmd.color("cyan", "receptor")
    cmd.color("salmon", "ligand")

    # Attempt to read grid center from param_dict. If missing, fallback = (0,0,0).
    try:
        cx = float(params_dict["grid_center_x"])
        cy = float(params_dict["grid_center_y"])
        cz = float(params_dict["grid_center_z"])
    except:
        cx, cy, cz = 0.0, 0.0, 0.0

    # Draw axes at the grid center
    add_axes(center=(cx, cy, cz), length=10.0, radius=0.2)

    # Draw bounding box if present
    try:
        sx = float(params_dict["grid_size_x"])
        sy = float(params_dict["grid_size_y"])
        sz = float(params_dict["grid_size_z"])
        draw_box((cx, cy, cz), (sx, sy, sz))
    except:
        pass

    # Camera
    cmd.center("ligand")
    cmd.zoom("ligand", 10)

    # Render
    cmd.ray(800, 600)
    cmd.png(png_output, width=800, height=600, dpi=72)
    logging.info(f"Saved docking image: {png_output}")

    cmd.save(pse_output)
    logging.info(f"Saved PyMOL session: {pse_output}")

    cmd.delete("all")


@logger_decorator
def draw_box(center, size):
    """
    Draws a simple bounding box (CGO) based on center=(cx, cy, cz) and size=(sx, sy, sz).
    """
    cx, cy, cz = center
    sx, sy, sz = size

    half_x = float(sx) / 2.0
    half_y = float(sy) / 2.0
    half_z = float(sz) / 2.0

    min_x = cx - half_x
    max_x = cx + half_x
    min_y = cy - half_y
    max_y = cy + half_y
    min_z = cz - half_z
    max_z = cz + half_z

    box = [
        LINEWIDTH, 2.0,
        COLOR, 1.0, 0.0, 1.0,
        BEGIN, LINES,

        # Bottom rectangle
        VERTEX, min_x, min_y, min_z, VERTEX, max_x, min_y, min_z,
        VERTEX, max_x, min_y, min_z, VERTEX, max_x, max_y, min_z,
        VERTEX, max_x, max_y, min_z, VERTEX, min_x, max_y, min_z,
        VERTEX, min_x, max_y, min_z, VERTEX, min_x, min_y, min_z,

        # Top rectangle
        VERTEX, min_x, min_y, max_z, VERTEX, max_x, min_y, max_z,
        VERTEX, max_x, min_y, max_z, VERTEX, max_x, max_y, max_z,
        VERTEX, max_x, max_y, max_z, VERTEX, min_x, max_y, max_z,
        VERTEX, min_x, max_y, max_z, VERTEX, min_x, min_y, max_z,

        # Vertical lines
        VERTEX, min_x, min_y, min_z, VERTEX, min_x, min_y, max_z,
        VERTEX, max_x, min_y, min_z, VERTEX, max_x, min_y, max_z,
        VERTEX, max_x, max_y, min_z, VERTEX, max_x, max_y, max_z,
        VERTEX, min_x, max_y, min_z, VERTEX, min_x, max_y, max_z,

        END
    ]
    cmd.load_cgo(box, "docking_box")


@logger_decorator
def generate_html_results(
    html_file,
    receptor_name,
    ligands_file,
    ligand_results,
    receptor_pdbqt
):
    """
    Generates an HTML file with docking results (no references to P2Rank).
    - 2D structure (SVG)
    - 3D docking image (PNG), clickable
    - Download links (.pdbqt, .pse)
    - Summaries sorted by best affinity
    - Link to text file with all docking poses
    """
    try:
        # Sort by best_affinity ascending
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
            hf.write('.docking-energy { background-color: #dfffe0; }\n')
            hf.write('</style>\n')
            hf.write('</head>\n')
            hf.write('<body>\n')

            receptor_pdbqt_basename = os.path.basename(receptor_pdbqt)
            header_text = (
                f"Docking results for receptor: "
                f'<span style="color: red;">{receptor_name}</span><br/>'
                f"Using ligand file: "
                f'<span style="color: navy;">{os.path.basename(ligands_file)}</span><br/>'
            )
            hf.write(f'<h2 style="text-align: center;">{header_text}</h2>\n')
            hf.write('<div style="text-align:center;">\n')
            hf.write(
                f'<a href="{receptor_pdbqt_basename}" download="receptor_structure.pdbqt" '
                f'type="application/octet-stream">Receptor structure (.PDBQT)</a>\n'
            )
            hf.write('</div><br/>\n')

            hf.write('<table>\n')
            hf.write(
                '<tr>'
                '<th>#</th>'
                '<th>Ligand Name</th>'
                '<th>2D Structure (SVG)</th>'
                '<th>Docking Preview Image</th>'
                '<th class="docking-energy">Affinity (kcal/mol)</th>'
                '<th>.pdbqt Output</th>'
                '<th>PyMOL Session</th>'
                '</tr>\n'
            )

            for idx, item in enumerate(ligand_results_sorted, start=1):
                lig_name  = item['name']
                aff       = item['affinity']
                aff_str   = f"{aff:.2f}" if aff is not None else 'N/A'

                # Paths relative to HTML location
                image_svg_rel  = os.path.relpath(item['image'],         os.path.dirname(html_file))
                dock_img_rel   = os.path.relpath(item['docking_image'], os.path.dirname(html_file))
                pdbqt_out_rel  = os.path.relpath(item['output_pdbqt'],  os.path.dirname(html_file))
                pse_rel        = os.path.relpath(item['pymol_session'], os.path.dirname(html_file))

                pdbqt_link_text = os.path.basename(item['output_pdbqt'])
                pse_link_text   = os.path.basename(item['pymol_session'])

                hf.write('<tr>\n')
                hf.write(f'<td>{idx}</td>\n')
                hf.write(f'<td>{lig_name}</td>\n')
                hf.write(f'<td><img src="{image_svg_rel}" alt="2D structure" style="width:300px;"/></td>\n')
                hf.write(
                    f'<td><a href="{dock_img_rel}" target="_blank">'
                    f'<img src="{dock_img_rel}" alt="Docking Image" style="max-width:150px; max-height:150px;"/></a></td>\n'
                )
                hf.write(f'<td class="docking-energy">{aff_str}</td>\n')
                hf.write(
                    f'<td><a href="{pdbqt_out_rel}" download="{pdbqt_link_text}" type="application/octet-stream">'
                    f'{pdbqt_link_text}</a></td>\n'
                )
                hf.write(
                    f'<td><a href="{pse_rel}" download="{pse_link_text}" type="application/octet-stream">'
                    f'{pse_link_text}</a></td>\n'
                )
                hf.write('</tr>\n')

            hf.write('</table>\n')
            hf.write('<br/>\n')

            # Link to the text file with all docking poses
            text_file_name = f"{receptor_name}_results.txt"
            hf.write('<p style="text-align:center;">\n')
            hf.write(
                f'<a href="{text_file_name}" target="_blank">'
                f'Full docking poses in a text file (all modes). CLICK</a>\n'
            )
            hf.write('</p>\n')

            hf.write('</body>\n')
            hf.write('</html>\n')

        logging.info(f"HTML report generated: {html_file}")
    except Exception as e:
        logging.error(f"Error in generating HTML results: {e}")


if __name__ == "__main__":
    main()