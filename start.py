import subprocess
import sys
import os

def run_command(command_list):
    try:
        subprocess.run(command_list, check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

# Install Python libraries
print("Installing Python libraries...")
run_command(['pip', 'install', '--upgrade', 'pip'])
run_command(['pip', 'install', 'biopandas', 'pubchempy', 'tqdm', 'matplotlib', 'scipy', 'rdkit-pypi', 'biopython'])

# Install system tools
print("Updating package list...")
run_command(['sudo', 'apt-get', 'update'])

print("Installing system packages...")
run_command(['sudo', 'apt-get', 'install', '-y', 'pymol', 'openbabel', 'wget', 'tar', 'default-jre'])

# Install PyMOL Python module if necessary
try:
    import pymol
except ImportError:
    print("PyMOL Python module not found. Installing...")
    run_command(['pip', 'install', 'pymol-open-source'])

# Download AutoDock Vina 1.2.5
print("Downloading AutoDock Vina 1.2.5...")
run_command(['wget', 'https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64'])

# Make the file executable
print("Making AutoDock Vina executable...")
run_command(['chmod', '+x', 'vina_1.2.5_linux_x86_64'])

# Move the file to /usr/local/bin
print("Moving AutoDock Vina to /usr/local/bin...")
run_command(['sudo', 'mv', 'vina_1.2.5_linux_x86_64', '/usr/local/bin/'])

# Check the installed version
print("Checking AutoDock Vina version...")
run_command(['vina_1.2.5_linux_x86_64', '--version'])

print("AutoDock Vina 1.2.5 installed successfully!")

# Download and extract p2rank
print("Downloading p2rank...")
run_command(['wget', 'https://github.com/rdk/p2rank/releases/download/2.4.2/p2rank_2.4.2.tar.gz'])

print("Extracting p2rank...")
run_command(['tar', '-xzf', 'p2rank_2.4.2.tar.gz'])

print("p2rank installed successfully!")
