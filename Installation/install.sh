#!/bin/bash

# Installation script for Ubuntu 22.04
# Downloads and installs Miniconda, creates the "auto_dock" environment, and runs the provided Python script

# Step 1: Download and install Miniconda
echo "Downloading Miniconda..."
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
echo "Installing Miniconda..."
bash miniconda.sh -b -p $HOME/miniconda
rm miniconda.sh

# Add Miniconda to PATH
export PATH="$HOME/miniconda/bin:$PATH"
echo "source \$HOME/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc
source $HOME/miniconda/etc/profile.d/conda.sh

# Step 2: Create a conda environment with Python 3.11
echo "Creating the conda environment: auto_dock..."
conda create -y -n auto_dock python=3.11

# Step 3: Activate the conda environment
echo "Activating the auto_dock environment..."
conda activate auto_dock

# Step 4: Run the provided Python script
echo "Running the provided Python script..."
python modules.py

echo "Script executed successfully!"
