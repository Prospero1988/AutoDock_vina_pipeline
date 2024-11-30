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

# Install Python libraries using Conda
print("Installing Python libraries with Conda...")
run_command(['conda', 'install', '-y', '-c', 'conda-forge',
             'biopython', 'biopandas', 'pubchempy', 'tqdm', 
             'matplotlib', 'scipy', 'rdkit', 'pdbfixer'])

# Ensure pip is installed in the current Conda environment
print("Ensuring pip is installed in Conda environment...")
run_command(['conda', 'install', '-y', '-c', 'conda-forge', 'pip'])
run_command(['conda', 'install', '-y', '-c', 'conda-forge', 'pymol-open-source'])

# Install remaining Python libraries with pip (if not available in Conda)
print("Installing additional Python libraries with pip...")
run_command(['pip', 'install', '--upgrade', 'pip'])
#run_command(['pip', 'install', 'pymol-open-source'])

# Install system tools
print("Updating package list...")
run_command(['sudo', 'apt-get', 'update'])

print("Installing system packages...")
run_command(['sudo', 'apt-get', 'install', '-y', 'openbabel', 'wget', 'tar', 'default-jre'])

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

# Download and install P2Rank
print("Downloading P2Rank...")
run_command(['wget', 'https://github.com/rdk/p2rank/releases/download/2.4.2/p2rank_2.4.2.tar.gz'])

print("Extracting P2Rank...")
run_command(['tar', '-xzf', 'p2rank_2.4.2.tar.gz'])

# Move P2Rank to /opt
print("Moving P2Rank to /opt/p2rank...")
if os.path.exists('/opt/p2rank'):
    print("Removing existing P2Rank installation...")
    run_command(['sudo', 'rm', '-rf', '/opt/p2rank'])
run_command(['sudo', 'mv', 'p2rank_2.4.2', '/opt/p2rank'])

# Adjust P2Rank directory structure
print("Adjusting P2Rank directory structure...")
run_command(['sudo', 'mv', '/opt/p2rank/*', '/opt/p2rank'])
run_command(['sudo', 'rmdir', '/opt/p2rank/p2rank_2.4.2'])

# Update permissions for P2Rank
print("Updating permissions for /opt/p2rank...")
run_command(['sudo', 'chmod', '-R', '755', '/opt/p2rank'])

# Update the prank script
print("Updating prank script...")
prank_script_path = '/opt/p2rank/prank'
if os.path.exists(prank_script_path):
    with open(prank_script_path, 'r') as f:
        prank_script = f.readlines()
    
    prank_script = [line.replace("$THIS_SCRIPT_DIR_REL_PATH", "/opt/p2rank") for line in prank_script]
    
    with open(prank_script_path, 'w') as f:
        f.writelines(prank_script)
    run_command(['sudo', 'chmod', '+x', prank_script_path])

# Create a symlink for prank
print("Creating symlink for prank...")
run_command(['sudo', 'ln', '-sf', '/opt/p2rank/prank', '/usr/local/bin/prank'])

print("P2Rank installed and configured successfully!")
