import subprocess
import sys
import os
import shutil

def run_command(command_list):
    try:
        subprocess.run(command_list, check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

def is_java_installed():
    try:
        subprocess.run(['java', '-version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        return True
    except FileNotFoundError:
        return False

# Check if Java is installed
if not is_java_installed():
    print("Java is not installed. Installing Java...")
    run_command(['sudo', 'apt-get', 'update'])
    run_command(['sudo', 'apt-get', 'install', '-y', 'default-jre'])
else:
    print("Java is already installed.")

# Install Python libraries using Conda
print("Installing Python libraries with Conda...")
run_command(['conda', 'install', '-y', '-c', 'conda-forge',
             'biopython', 'biopandas', 'pubchempy', 'tqdm', 
             'matplotlib', 'scipy', 'rdkit', 'pdbfixer', 'pymol-open-source'])

# Ensure pip is installed in the current Conda environment
print("Ensuring pip is installed in Conda environment...")
run_command(['conda', 'install', '-y', '-c', 'conda-forge', 'pip'])
#run_command(['conda', 'install', '-y', '-c', 'conda-forge', 'pymol-open-source'])

# Install remaining Python libraries with pip (if not available in Conda)
print("Installing additional Python libraries with pip...")
run_command(['pip', 'install', '--upgrade', 'pip'])

# Install system tools
print("Updating package list...")
run_command(['sudo', 'apt-get', 'update'])

print("Installing system packages...")
run_command(['sudo', 'apt-get', 'install', '-y', 'openbabel', 'wget', 'tar'])
