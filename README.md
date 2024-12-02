
# Docking System

This repository provides an automated docking solution for ligands and receptor proteins using AutoDock Vina and P2Rank. It supports high-throughput docking workflows and integrates seamlessly with SLURM or can be run locally.

## Technologies Used
- **Python 3.11**: Core scripting language.
- **AutoDock Vina v1.2.5**: Molecular docking engine.
- **P2Rank v2.4.2**: Binding site prediction.
- **Biopython, RDKit, Open Babel, PyMOL**: Molecular handling, visualization, and preparation tools.
- **SLURM**: Workload manager for distributed computing (optional).

## Requirements
- **Ubuntu 22.04**
- **Miniconda** (Installed via provided scripts)
- **SLURM** (Optional for distributed execution)

### Python Libraries
- `biopython`, `biopandas`, `pubchempy`, `tqdm`, `matplotlib`, `scipy`, `rdkit`, `pdbfixer`, `pymol-open-source`

### System Tools
- `openbabel`, `wget`, `tar`
- **Java Runtime Environment (JRE)**

## Installation

### Full Installation (Fresh System)
1. Clone the repository:
    ```bash
    git clone https://github.com/your-repository/docking-system.git
    cd docking-system
    ```
2. Run the installation script:
    ```bash
    chmod +x install.sh
    bash install.sh
    ```

### Minimal Installation (Configured System)
For environments where most dependencies are already configured:
```bash
chmod +x mini_install.sh
bash mini_install.sh
```

### Additional Configuration
Ensure the following tools are available in their respective paths:
- **AutoDock Vina**: `/usr/local/bin/vina_1.2.5_linux_x86_64`
- **P2Rank**: `/usr/local/bin/prank`

## Usage

### SLURM Execution
1. Prepare input files:
   - Place receptor PDB IDs in a CSV file under `./receptors`.
   - Place ligand structures in SDF format under `./ligands`.

2. Submit the job via SLURM:
    ```bash
    sbatch start_docking.sh
    ```

### Local Execution
1. Activate the conda environment:
    ```bash
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate auto_dock
    ```
2. Run the Python script:
    ```bash
    python3 init_docking.py --pdb_ids receptors.csv --ligands ligand_file.sdf
    ```

## SLURM Configuration
The repository includes a sample SLURM script (`start_docking.sh`) optimized for the docking pipeline. Key configurations include:
- Single task allocation (`#SBATCH --ntasks=1`).
- Infinite runtime (`#SBATCH --time=INFINITE`).

## Input Parameters
- `--pdb_ids`: CSV file with receptor PDB codes.
- `--ligands`: SDF file containing ligands.
- `--tol`: Docking box tolerance (Å, default: 0).
- `--pckt`: Pocket number from P2Rank predictions (default: 1).
- `--exhaust`: Docking thoroughness (default: 20).
- `--energy_range`: Energy range for docking poses (default: 2 kcal/mol).

## Outputs
1. Results organized by receptor:
   - `receptor_name_results.txt`: Detailed docking results.
   - `ligand_name.pdbqt`: Prepared ligand.
   - `ligand_name.svg`: 2D ligand structure.

2. HTML Report:
   - Summarized docking results.
   - Interactive visualization links.

## Notes
- The system works best with SLURM for distributed execution but can run locally.
- Ensure all dependencies are correctly installed and configured.
- Follow the user manual (`User_Guide_Docking_System_ENG.html`) for detailed steps.

For more details, refer to the [Installation Guide](Installation_Guide_ENG.html).
