# AutoDock Vina Pipeline with P2RANK for HTS Docking with WEB Interface

## Table of Contents

- [Description](#description)
- [UPDATES](#updates)
- [Web Interface](#web-interface)
- [Output Example](#output-example)
- [Results Screenshots](#results-screenshots)
- [Technologies Used and Requirements](#technologies-used)
- [Installation](#installation)
  - [1. Create a User Account for Docking](#1-create-a-user-account-for-docking)
  - [2. Clone the Repository](#2-clone-the-repository)
  - [3. Install Dependencies and Configure](#3-install-dependencies-and-configure)
  - [4. Install and Configure SLURM for Task Management](#4-install-and-configure-slurm-for-task-management)
  - [5. Configure the Server as a System Service](#5-configure-the-server-as-a-system-service)
  - [6. NGINX Reverse Proxy Setup for Streamlit](#6-nginx-reverse-proxy-setup-for-streamlit)
- [Running the Application](#running-the-application)
  - [Accessing Locally](#accessing-locally)
  - [Accessing on LAN](#accessing-on-lan)
  - [Functionalities Overview](#functionalities-overview)
  - [Instructions for Use](#instructions-for-use)
- [Docking Work Flow](#docking-work-flow)
  - [Input Parsing](#input-parsing)
  - [Receptor Preparation](#receptor-preparation)
  - [Binding Site Prediction](#binding-site-prediction)
  - [Ligand Preparation](#ligand-preparation)
  - [Docking Execution](#docking-execution)
- [Docking Parameters explained](#docking-parameters-explained)
  - [Tolerance parameters](#tolerance-parameters)
  - [Auto Dock Vina parameters](#auto-dock-vina-parameters)
- [Visualization and Results Generation](#visualization-and-results-generation)
  - [Visualizations](#visualizations)
  - [HTML and CSV Reports](#html-and-csv-reports)
- [Output File Organization](#output-file-organization)
- [Acknowledgments](#acknowledgments)

## Description

### Program Description

This repository provides an automated docking solution for ligands and receptor proteins using **AutoDock Vina** and **P2Rank**. It supports high-throughput docking workflows and integrates seamlessly with **SLURM**, a workload manager for distributed computing and queue management. Additionally, the entire functionality is built into a web interface using the **Streamlit** framework, allowing convenient and intuitive operation from a web browser.

The program is built on open-source libraries and solutions. It implements a user account system, facilitating easy project management, handling of generated data, and seamless navigation. The installation process is comprehensive, enabling even beginner users to utilize this tool effectively. The system is developed and tested on **Ubuntu 22.04** as a compute hosting server. Access via a web browser is unrestricted by the system, allowing the program to run locally on an Ubuntu machine or be configured as a server on a local network, accessible from any computer within the same LAN.

### Notes

I am still working on this project and adding new features. The uploaded repository is fully functional. If you have any suggestions then feel free to contact me, or open a discussion or add a post.

### Problems and support

If you encounter any problems with the installation or operation of the docking program, do not hesitate to contact me. I will do my best to assist you.

## UPDATES

- Some PDB codes in the PDB database do not correspond to actual .pdb files, because only mmCIF files are available. The program, in case it cannot download a .pdb file, tries to download a .cif file. After successfully downloading it, it converts it to .pdb and passes it on for further processing and docking.
- Automatic selection of a chain with a receptor did not work due to the lack of uniformity in chain naming. Now, when PDB codes are given, the structures are retrieved, a list of chains is loaded, and the user has to select for each receptor the appropriate chain containing the receptor, or docking site, from a drop-down list. Alternatively, an input can be prepared in the form of a csv file that contains PDB codes and chain IDs. After loading the CSV file, a list of receptors with selected chains is shown, which can be modified, or accepted and passed on for further calculations.

## Web Interface

The interface is built on the Sreamlit framework. After installation and configuration, the whole thing functions perfectly from a web browser.

</br>
<div align="center">
    <img src="img/interface.png" alt="Presentation of web interface" title="Presentation of web interface">
</div>

## Output Example

In this repository, there is a downloadable file [output_example.zip](https://github.com/Prospero1988/AutoDock_vina_pipeline/raw/refs/heads/main/output_example.zip) which contains an archive with a zipped sample output for the docking program.

## Results Screenshots
</br>
<div align="center">
    <img src="img/results.png" alt="Docking results visualization" title="Visualization of docking results">
</div>
</br>
<div align="center">
    <img src="img/grid.png" alt="Docking results visualization" title="Visualization of docking results">
</div>
</br>
<div align="center">
    <img src="img/csv.png" alt="Docking results visualization" title="Visualization of docking results">
</div>
</br>
<hr>

## Technologies Used

- **Python 3.11**: Core scripting language.
- **AutoDock Vina v1.2.5**: Molecular docking engine.
- **P2Rank v2.4.2**: Binding site prediction.
- **Biopython**, **RDKit**, **Open Babel**, **open-PyMOL**: Molecular handling, visualization, and preparation tools.
- **SLURM**: Workload manager for distributed computing and queue management.
- **Streamlit 1.40.2**: Frontend web interface.

## Requirements

- **Operating System**: Ubuntu 22.04 (or other compatible Debian distributions). For advanced users, any Linux distribution can be used, but library and installation package adjustments may be necessary.
  
### Python Libraries

- `biopython`
- `biopandas`
- `pandas`
- `pubchempy`
- `tqdm`
- `matplotlib`
- `scipy`
- `rdkit`
- `pdbfixer`
- `pymol-open-source`
- `streamlit`
- `bcrypt`

**Java Runtime Environment (JRE)**

## Installation

The installation process is divided into several main stages. The program is configured to operate under a specific user account and name. If you wish to modify this, locate all instances of `docking_machine` in the `dock_GUI.py` file and replace them with your desired username. The docking server will be set up on this account.

### 1. Create a User Account for Docking

Create a user account named `docking_machine` and assign the necessary administrative permissions to install required packages and libraries. These permissions can be revoked after installation.

```bash
sudo adduser docking_machine
sudo usermod -aG sudo docking_machine
```

Log into the new user account:

```bash
su - docking_machine
```

*Alternatively, relog into this account by opening a new terminal session if you are working via SSH.*

### 2. Clone the Repository

Clone the repository and move it to the preferred directory, e.g., `/home/docking_machine/dock`:

```bash
git clone https://github.com/Prospero1988/AutoDock_vina_pipeline.git
mkdir ~/dock
mv ~/AutoDock_vina_pipeline/* ~/dock/
mv ~/AutoDock_vina_pipeline/.* ~/dock/ 2>/dev/null
rmdir ~/AutoDock_vina_pipeline
```

*Alternatively, you can perform these steps manually using your operating system's graphical interface.*

Set access permissions for the `dock` directory:

```bash
chmod -R 755 /home/docking_machine/dock
sudo chown -R docking_machine:docking_machine /home/docking_machine/dock
```

### 3. Install Dependencies and Configure

Navigate to the installation directory and run the installation script. **Do not** start this script with `sudo` to ensure that all packages and libraries are installed under the `docking_machine` user account, not the root account.

```bash
cd /home/docking_machine/dock/installation
chmod +x install.sh
bash install.sh
```

At this stage, you can already work with the program via the command line by running the `init_docking.py` Python script. However, this method is not very convenient and does not support the workload manager or the graphical interface/server setup.

### 4. Install and Configure SLURM for Task Management

Ensure the system is up to date:

```bash
sudo apt update && sudo apt upgrade -y
```

Install the required dependencies:

```bash
sudo apt install -y munge libmunge-dev libmunge2 build-essential slurm-wlm slurm-client
```

Create an authentication key for Munge:

```bash
sudo /usr/sbin/create-munge-key
```

If it's not working you have to create munge-key manually. Log into root account:

```bash
su - root
```

Create munge key and add permissions:

```bash
sudo dd if=/dev/urandom bs=1 count=1024 > /etc/munge/munge.key
sudo chown munge:munge /etc/munge/munge.key
sudo chmod 400 /etc/munge/munge.key
```

Check if key was created:

```bash
ls -l /etc/munge/munge.key
```

You should see something like this:
```bash
-r-------- 1 munge munge 1024 <date of creation> /etc/munge/munge.key
```

Set appropriate permissions:

```bash
sudo chown -R munge: /etc/munge /var/lib/munge /var/log/munge
sudo chmod 700 /etc/munge /var/lib/munge /var/log/munge
```

Start and enable the Munge service:

```bash
sudo systemctl enable munge
sudo systemctl start munge
```

Verify Munge is working correctly:

```bash
munge -n | unmunge
```

*Expected output: Success (0).*

Add a dedicated user for SLURM:

```bash
sudo useradd -r -m -d /var/lib/slurm -s /bin/false slurm
```

Create necessary directories:

```bash
sudo mkdir -p /var/spool/slurmd /var/log/slurm
```

Set appropriate permissions:

```bash
sudo chown -R slurm: /var/spool/slurmd /var/log/slurm
sudo chmod -R 755 /var/spool/slurmd /var/log/slurm
```

Check the hostname:

```bash
hostname
```

Edit the SLURM configuration file:

```bash
sudo nano /etc/slurm/slurm.conf
```

Add the following minimal configuration, replacing `<YOUR_CLUSTER_NAME>` with your preferred cluster name and `<YOUR_HOSTNAME>` with the hostname obtained earlier:

```conf
# Basic Configuration
ClusterName=<YOUR_CLUSTER_NAME>
ControlMachine=<YOUR_HOSTNAME>

# Ports and Authentication
SlurmctldPort=6817
SlurmdPort=6818
AuthType=auth/munge

# Logging
SlurmdLogFile=/var/log/slurm/slurmd.log
SlurmctldLogFile=/var/log/slurm/slurmctld.log

# Resource Management
SlurmUser=slurm
StateSaveLocation=/var/spool/slurmd
SlurmdSpoolDir=/var/spool/slurmd
ProctrackType=proctrack/pgid
TaskPlugin=task/none
SchedulerType=sched/backfill

# Node Configuration
NodeName=localhost CPUs=16 RealMemory=64000 State=UNKNOWN

# Partition Configuration
PartitionName=main Nodes=ALL Default=YES MaxTime=INFINITE State=UP
```

Save the file and set appropriate permissions:

```bash
sudo chown slurm: /etc/slurm/slurm.conf
sudo chmod 644 /etc/slurm/slurm.conf
```

Start and enable SLURM services:

```bash
sudo systemctl enable slurmctld
sudo systemctl enable slurmd
sudo systemctl start slurmctld
sudo systemctl start slurmd
```

Check the status of SLURM services:

```bash
sudo systemctl status slurmctld
sudo systemctl status slurmd
```

*Expected output: Active: active (running).*

Ensure all services start automatically on boot:

```bash
sudo systemctl enable munge
sudo systemctl enable slurmctld
sudo systemctl enable slurmd
```

### 5. Configure the Server as a System Service

To ensure the server is always available after a system restart, configure it as a system service using **systemd**.

Create a `streamlit_docking.service` file in the `/etc/systemd/system` directory:

```bash
sudo nano /etc/systemd/system/streamlit_docking.service
```

Add the following configuration, replacing paths as necessary:

```ini
[Unit]
Description=Streamlit Docking GUI Service
After=network.target

[Service]
User=docking_machine
Group=docking_machine
WorkingDirectory=/home/docking_machine/dock
Environment="PATH=/home/docking_machine/miniconda/envs/auto_dock/bin"
ExecStart=/home/docking_machine/miniconda/envs/auto_dock/bin/streamlit run dock_GUI.py
Restart=always
RestartSec=5

[Install]
WantedBy=multi-user.target
```

Set the correct permissions and ownership:

```bash
sudo chown root:root /etc/systemd/system/streamlit_docking.service
sudo chmod 644 /etc/systemd/system/streamlit_docking.service
```

Enable and start the service:

```bash
sudo systemctl daemon-reload
sudo systemctl enable streamlit_docking.service
sudo systemctl start streamlit_docking.service
```

Verify the service status:

```bash
sudo systemctl status streamlit_docking.service
```

*Expected output:*
```
● streamlit_docking.service - Streamlit Docking GUI Service
     Loaded: loaded (/etc/systemd/system/streamlit_docking.service; enabled; vendor preset: enabled)
     Active: active (running) since ...
```

Open the necessary ports in the firewall:

```bash
sudo ufw allow 8001,8501
sudo ufw status
```

Test the setup after a system restart:

```bash
sudo reboot
```

After rebooting, verify the service is running:

```bash
sudo systemctl status streamlit_docking.service
```

Access the application via a web browser at:

```
http://<YOUR_COMPUTER_IP>:8501
```

### 6. NGINX Reverse Proxy Setup for Streamlit

To access the Streamlit application using the address `http://<YOUR_LAN_IP>` without specifying a port, you need to configure a reverse proxy using a web server like **NGINX**. The reverse proxy will handle requests on the default HTTP port (80) and forward them to Streamlit on port 8501.

#### Instructions for NGINX

##### 1. Install NGINX (if not already installed)

```bash
sudo apt update
sudo apt install nginx
```

##### 2. Configure NGINX as a Reverse Proxy

Create a configuration file for your Streamlit application:

```bash
sudo nano /etc/nginx/sites-available/streamlit
```

Insert the following configuration, replacing `<YOUR_LAN_IP>` with your actual LAN IP address:

```nginx
server {
    listen 80;
    server_name <YOUR_LAN_IP>;

    location / {
        proxy_pass http://127.0.0.1:8501;  # Forward requests to Streamlit
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;

        # Handle WebSocket connections (required for Streamlit)
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
    }
}
```

##### 3. Enable the Configuration

Create a symbolic link to the `sites-enabled` directory:

```bash
sudo ln -s /etc/nginx/sites-available/streamlit /etc/nginx/sites-enabled/
```

##### 4. Remove Default Configuration (Optional)

If the default NGINX configuration is not needed, remove it:

```bash
sudo rm /etc/nginx/sites-enabled/default
```

##### 5. Test Configuration and Restart NGINX

Test the NGINX configuration to ensure there are no errors:

```bash
sudo nginx -t
```

If everything is correct, restart NGINX:

```bash
sudo systemctl restart nginx
```
##### 6. Open Firewall Ports for NGINX

To allow HTTP and HTTPS traffic through the firewall, use the following command:

```bash
sudo ufw allow 'Nginx Full'
```

Verify the firewall status:

```bash
sudo ufw status
```

Expected output:

```bash
To                         Action      From
--                         ------      ----
Nginx Full                 ALLOW       Anywhere
Nginx Full (v6)            ALLOW       Anywhere (v6)
```

##### 7. Test Access

Open your browser and navigate to:

```
http://<YOUR_LAN_IP>
```

#### NGINX as a System Service

NGINX is designed as a system service that automatically starts during system boot. By default, after installation, NGINX is configured to start automatically. However, you can verify or manually enable this option.

##### Check NGINX Status

Run the following command to ensure NGINX is enabled and starts on boot:

```bash
sudo systemctl is-enabled nginx
```

The expected output is:

```
enabled
```

##### Enable NGINX at Startup (If Not Enabled)

If the previous command returns `disabled`, you can enable NGINX:

```bash
sudo systemctl enable nginx
```

##### Verify NGINX is Running

Check if NGINX is active:

```bash
sudo systemctl status nginx
```

If everything is fine, you should see a message similar to:

```
● nginx.service - A high performance web server and a reverse proxy server
     Loaded: loaded (/lib/systemd/system/nginx.service; enabled; vendor preset: enabled)
     Active: active (running) since ...
```

#### Summary

By configuring NGINX as a service and enabling it, it will run automatically at every system boot without needing manual intervention. If your Streamlit server is properly running through NGINX, accessing the application via `http://<YOUR_LAN_IP>` will consistently work after each restart. 



## Running the Application

### Accessing Locally

To access the application on the local computer where the server is running, open your web browser and navigate to:

```
http://localhost:8501
```

*Alternatively, you can use:*

```
http://127.0.0.1:8501
```

### Accessing on LAN

To access the application from other computers within your Local Area Network (LAN):

1. **Ensure NGINX is Properly Configured**: As outlined in the [NGINX Reverse Proxy Setup for Streamlit](#6-nginx-reverse-proxy-setup-for-streamlit) section, NGINX should be set up to forward requests from port 80 to Streamlit's port 8501.

2. **Open Firewall Ports**: Ensure that ports 80 (HTTP) and 443 (HTTPS) are open in your firewall to allow incoming traffic. This was covered in step 6 of the NGINX setup.

3. **Determine Server's LAN IP Address**: Find the LAN IP address of the server hosting the application by running:

    ```bash
    hostname -I
    ```

    Suppose the LAN IP is `192.168.1.100`.

4. **Access from Other Devices**: On any device within the same LAN, open a web browser and navigate to:

    ```
    http://192.168.1.100
    ```

    Replace `192.168.1.100` with your server's actual LAN IP address.

### Functionalities Overview

The Docking Program provides the following functionalities:

- **User Authentication**: Secure login and registration system with password hashing using bcrypt.

- **Project Management**: Create, manage, and delete docking projects, each containing receptor proteins and ligand files.

- **Docking Submission**: Configure docking parameters and submit docking jobs to the SLURM workload manager for high-throughput molecular docking.

- **Queue Monitoring**: View the current SLURM job queue and monitor the status of submitted jobs.

- **Results Visualization and Download**: Access interactive results generated by the docking process and download them as CSV files.

- **PyMOL Installation Guide**: Access a guide for installing Open-PyMOL on Windows for molecular visualization.

### Instructions for Use

1. **Login/Register**:
   - Open the application in your web browser using either `http://localhost:8501` or `http://<YOUR_LAN_IP>`.
   - If you already have an account, enter your username and password and click "Login".
   - To create a new account, click "Add New User", enter your desired username and password, and click "Register".

2. **Main Menu**:
   - Upon successful login, you will be greeted with the main menu, where you can select from the following modules:
     - **DOCKING**: Configure and submit docking jobs.
     - **QUEUE**: View and manage SLURM job queue.
     - **SHOW RESULTS**: View interactive results of docking simulations.
     - **DOWNLOAD RESULTS**: Download docking results.
     - **DELETE RESULTS**: Delete existing projects and results.
     - **PyMOL Installation GUIDE**: Access the guide for installing Open-PyMOL on Windows.
     - **LOG OUT**: Log out of the application.

3. **DOCKING Module**:
   - **Project Setup**: Enter a project name to create a new docking project.
   - **Enter PDB Codes**: Provide PDB codes for receptor proteins, either by typing them in or uploading a CSV file.
   - **Upload Ligand Files**: Upload ligand files in `.mol2` or `.SDF` formats.
   - **Docking Parameters**: Configure docking parameters as needed.
   - **Project Summary**: Review your project details before submitting the docking job.
   - **Start Docking**: Submit the docking job to SLURM for processing.

4. **QUEUE Module**:
   - View the current SLURM job queue, including job IDs, users, job names, states, time used, and start times.
   - Refresh the queue status as needed.

5. **SHOW RESULTS Module**:
   - Select a project and receptor to view interactive docking results.
   - Access interactive results in a new browser tab.
   - Download results as CSV files.

6. **DOWNLOAD RESULTS Module**:
   - Select projects to download their results as a ZIP archive.

7. **DELETE RESULTS Module**:
   - Select projects to delete their associated data and results.

8. **PyMOL Installation GUIDE**:
   - Access a detailed guide for installing Open-PyMOL on Windows to visualize molecular structures.

9. **LOG OUT**:
   - Log out of your account securely.

### Executable Python Script

The Streamlit application is defined in `dock_GUI.py`. Ensure that this script is properly installed and configured during the installation process (as described in the Installation section). The Streamlit server is managed as a system service, so it starts automatically on system boot.

You can also manually run the application by navigating to the installation directory and executing:

```bash
streamlit run dock_GUI.py
```

However, it's recommended to use the system service configuration for seamless operation.

## Docking Work Flow

### Input Parsing

The script accepts the following arguments:

- `--pdb_ids`: A CSV file located in the `./receptors` directory, containing the PDB IDs of receptor proteins. Each ID corresponds to a unique protein structure available in the Protein Data Bank (PDB).
- `--ligands`: A ligand file located in the `./ligands` directory. Supported formats include SDF and MOL2 files, allowing flexibility in ligand input.

**Tolerance Parameters:**

- `--tol_x`: Tolerance in Ångströms to expand the docking pocket dimension along the X-axis beyond those defined by P2Rank (default: 0.0).
- `--tol_y`: Tolerance in Ångströms to expand the docking pocket dimension along the Y-axis beyond those defined by P2Rank (default: 0.0).
- `--tol_z`: Tolerance in Ångströms to expand the docking pocket dimension along the Z-axis beyond those defined by P2Rank (default: 0.0).

**Offset Parameters:**

- `--offset_x`: Offset in Ångströms to shift the center of the docking grid box along the X-axis (optional, default: 0.0).
- `--offset_y`: Offset in Ångströms to shift the center of the docking grid box along the Y-axis (optional, default: 0.0).
- `--offset_z`: Offset in Ångströms to shift the center of the docking grid box along the Z-axis (optional, default: 0.0).
- `--pckt`: Pocket number to use from P2Rank predictions (default: 1).
- `--exhaust`: Specifies how thorough the search should be for the best binding poses. Higher values increase precision but require more computation time (default: 16).
- `--energy_range`: Determines the range of energy scores (in kcal/mol) for poses to be considered (default: 4).

**Automatic Ligand Naming:**

In cases where ligands in the input files lack explicit names, the script assigns them generic names in the format `ligand_001`, `ligand_002`, etc., ensuring consistent and organized output.

### Receptor Preparation

**Download Receptor Structures:**

For each PDB ID listed in the CSV file, the script downloads the corresponding protein structure from the Protein Data Bank (PDB). The downloaded file is saved as `<PDB_ID>_dirty.pdb` in a newly created folder named after the receptor (e.g., `./8W88/`).

**Fixing the Receptor:**

Using PDBFixer, the script:

- Retains only the chain with the maximum number of residues.
- Removes heteroatoms and water molecules.
- Adds missing residues, atoms, and hydrogens based on a physiological pH of 7.4.

The fixed structure is saved as `<PDB_ID>_fixed.pdb`.

**Receptor Conversion:**

The fixed PDB structure is converted to the `.pdbqt` format required by AutoDock Vina. The converted file is saved as `<PDB_ID>.pdbqt`.

### Binding Site Prediction

- The script utilizes **P2Rank** to predict potential binding sites (pockets) on the receptor.
  - The predictions are saved in a folder named `01_p2rank_output` within the receptor's directory.
  - A CSV file (`<PDB_ID>_predictions.csv`) lists each pocket's coordinates, size, and scores.
- The selected pocket (based on the `--pckt` argument) is used to define the docking box dimensions. This includes the center coordinates (`center_x`, `center_y`, `center_z`) and sizes (`size_x`, `size_y`, `size_z`) with optional tolerances (`--tol_x`, `--tol_y`, `--tol_z`) for each axis. Additionally, optional offsets (`--offset_x`, `--offset_y`, `--offset_z`) allow for independent shifting of the docking grid center along each axis.

### Ligand Preparation

- For each ligand in the provided SDF or MOL2 file:
  - **Format Handling:**
    - The script automatically detects the file format (SDF or MOL2) and processes accordingly.
  - **Conversion and Processing:**
    - The ligand is converted to `.pdb` format using RDKit.
    - Hydrogen atoms are added, and a 3D conformer is generated for the ligand.
    - The `.pdb` file is converted to `.pdbqt` format required for docking using Open Babel.
- **Output Organization:**
  - The prepared ligand files are stored in the `02_ligands_results` subdirectory within the receptor's folder (e.g., `./8W88/02_ligands_results/`).

### Docking Execution

- The script runs **AutoDock Vina** for each receptor-ligand pair:
  - The docking box is defined using P2Rank predictions, with optional tolerances and offsets.
  - Parameters such as `--exhaust` (exhaustiveness) and `--energy_range` control the thoroughness and energy tolerance for pose scoring.
  - Docking results are saved in `.pdbqt` format, and key details (e.g., binding affinities) are extracted from the output.
- **Post-Docking File Management:**
  - All `.pdbqt` files for ligands after docking are collectively copied into the `03_ligands_PDBQT` folder, facilitating easy access without navigating through individual folders.

## Docking Parameters Explained

You can run the docking without changing/defining these parameters. Default values will be used. To view the docking box, it's best to open any resulting `.pse` file in PyMOL. This file contains a drawn grid box along with the coordinate axes.

### Tolerance Parameters

- `--tol_x`, `--tol_y`, `--tol_z`: Docking box tolerances in Ångströms to expand the docking pocket dimensions along the X, Y, and Z axes respectively (default: 0.0 for each). Define how much to expand the docking pocket along each respective axis beyond the predictions made by P2Rank. Negative values indicate contraction.
- `--offset_x`, `--offset_y`, `--offset_z`: Offsets in Ångströms to shift the center of the docking grid box along the X, Y, and Z axes respectively (default: 0.0 for each).

### Auto Dock Vina Parameters

- `--pckt`: Pocket number from P2Rank predictions (default: 1).
- `--exhaust`: Docking thoroughness (default: 16).
- `--energy_range`: Energy range for docking poses (default: 4 kcal/mol).

## Visualization and Results Generation

### Visualizations

PyMOL is used to generate visualizations of the best-docked ligand poses superimposed on the receptor structure. The visualizations include the docking grid (box) and XYZ axes for better spatial orientation. Both high-resolution images (`.png`) and PyMOL session files (`.pse`) are saved for each docking result. The PyMOL session files include the docking grid and XYZ axes, allowing users to explore the docking results interactively within PyMOL.

### HTML and CSV Reports

#### HTML Report

The script creates an interactive HTML report for each receptor, summarizing:

- Key docking metrics (binding energies, pocket scores).
- Links to output files (e.g., `.pdbqt` and `.txt`).
- 2D and 3D visualizations of ligand-receptor complexes.
- **New**: An additional column with links to the PyMOL session files (`.pse`), enabling users to open and manipulate the docking results directly in PyMOL.

#### CSV Summary

In addition to the HTML report, a CSV summary file is generated containing:

- **Name**: Ligand name.
- **Affinity**: Binding affinity values.
- **SMILES**: Simplified molecular-input line-entry system representations of ligands.

This CSV file provides a convenient overview of docking results for further analysis.

## Output File Organization

In each project directory, there are subfolders corresponding to each receptor. Within these receptor-specific subfolders, you can find the results from the docking simulations. The organization of these subfolders is as follows:

Each receptor has its dedicated directory containing:

**Processed Structures:**

- `<PDB_ID>_dirty.pdb`: Raw receptor structure.
- `<PDB_ID>_fixed.pdb`: Cleaned receptor structure.
- `<PDB_ID>.pdbqt`: Receptor ready for docking.

**Docking Results:**

- `02_ligands_results/`:
  - `<ligand_name>.pdbqt`: Prepared ligand.
  - `<ligand_name>.svg`: 2D ligand structure images.
  - `<PDB_ID>_<ligand_name>_docking.pse`: PyMOL session files including the docking grid and XYZ axes.
- `03_ligands_PDBQT/`:
  - All docked ligand `.pdbqt` files copied here for easy access.

**Visualizations:**

- `<PDB_ID>_<ligand_name>_docking.png`: High-resolution 3D visualizations of docked complexes.

**Reports:**

- `<PDB_ID>_results.html`: Interactive HTML report summarizing docking results, including links to PyMOL session files.
- `<PDB_ID>_results_in_CSV.csv`: CSV file containing ligand name, affinity, and SMILES.

**P2Rank Predictions:**

- `01_p2rank_output/<PDB_ID>_predictions.csv`: Binding site information.

## Acknowledgments

- [AutoDock Vina](http://vina.scripps.edu/)
- [P2Rank](https://github.com/rdk/p2rank)
- [RDKit](https://www.rdkit.org/)
- [Open Babel](https://openbabel.org/)
- [PyMOL](https://pymol.org/)
- [Streamlit](https://streamlit.io/)
- [NGINX](https://nginx.org/en/)
- [SLURM](https://slurm.schedmd.com/overview.html)
---
