#!/bin/bash
#SBATCH --job-name=docking_job        
#SBATCH --output=docking_output.log
#SBATCH --error=docking_error.log
#SBATCH --ntasks=1                    
#SBATCH --time=INFINITE               
#SBATCH --partition=main              

# Activation of the conda environment
source ~/miniconda/etc/profile.d/conda.sh
conda activate auto_dock

# Running a Python script
python3 init_docking.py --pdb_ids 'receptors.csv' --ligands 'ligands.sdf' --tol 0 --pckt 1 --exhaust 16 --energy_range 3
