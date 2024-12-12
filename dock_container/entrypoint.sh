#!/bin/bash

# Start SLURM
slurmctld -c -f /etc/slurm/slurm.conf &
slurmd -c -f /etc/slurm/slurm.conf &

# Activate the Conda environment
source /root/miniconda/etc/profile.d/conda.sh
conda activate auto_dock

# Start Streamlit
streamlit run /root/dock/dock_GUI_container.py