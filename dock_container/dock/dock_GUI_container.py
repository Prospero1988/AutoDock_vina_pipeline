import streamlit as st
import streamlit.components.v1 as components
import os
import shutil
import re
import pandas as pd
import subprocess
import zipfile
import bcrypt
import threading
from http.server import SimpleHTTPRequestHandler, HTTPServer
from functools import partial
import requests
from Bio.PDB import PDBParser
import gzip


# Configure the Streamlit application
st.set_page_config(page_title="Docking Program", layout="centered")

# Define paths and constants
OBABEL_PATH = "/usr/bin/obabel"
HTTP_SERVER_PORT = 8001
RESULTS_DIR = '/root/dock/results'
os.makedirs(RESULTS_DIR, exist_ok=True)


@st.cache_resource
def start_shared_http_server():
    """
    Start a shared HTTP server that serves the 'static' directory for all users.

    Returns:
        HTTPServer: The running HTTPServer instance.
    """
    handler = partial(SimpleHTTPRequestHandler, directory='static')
    httpd = HTTPServer(('0.0.0.0', HTTP_SERVER_PORT), handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    return httpd


def validate_project_name(name):
    """
    Validate and sanitize a given project name.

    Args:
        name (str): The original project name.

    Returns:
        str or None: Returns a sanitized project name if valid, otherwise None.
    """
    name = name.replace(' ', '_')
    if re.match(r'^[A-Za-z0-9_]+$', name):
        return name
    else:
        return None


def hash_password(password):
    """
    Hash a password using bcrypt.

    Args:
        password (str): The plaintext password.

    Returns:
        str: The hashed password as a UTF-8 string.
    """
    salt = bcrypt.gensalt()
    hashed = bcrypt.hashpw(password.encode('utf-8'), salt)
    return hashed.decode('utf-8')


def check_credentials(username, password):
    """
    Check user credentials against stored hashes.

    Args:
        username (str): The username to check.
        password (str): The plaintext password to verify.

    Returns:
        bool: True if credentials match, otherwise False.
    """
    try:
        credentials = pd.read_csv(f'{RESULTS_DIR}/passwords.pw')
        user_row = credentials[credentials['user'] == username]
        if not user_row.empty:
            stored_hash = user_row.iloc[0]['password']
            if bcrypt.checkpw(password.encode('utf-8'), stored_hash.encode('utf-8')):
                return True
    except Exception as e:
        print(e)
        return False
    return False


def username_exists(username):
    """
    Check if a given username already exists in the system.

    Args:
        username (str): The username to check.

    Returns:
        bool: True if the username exists, otherwise False.
    """
    try:
        credentials = pd.read_csv(f'{RESULTS_DIR}/passwords.pw')
        return username in credentials['user'].values
    except Exception:
        return False


def add_new_user(username, hashed_password):
    """
    Add a new user with a hashed password to the system.

    Args:
        username (str): The new username.
        hashed_password (str): The hashed password.
    """
    try:
        if os.path.exists(f'{RESULTS_DIR}/passwords.pw'):
            credentials = pd.read_csv(f'{RESULTS_DIR}/passwords.pw')
            new_user = pd.DataFrame({'user': [username], 'password': [hashed_password]})
            credentials = pd.concat([credentials, new_user], ignore_index=True)
        else:
            credentials = pd.DataFrame({'user': [username], 'password': [hashed_password]})
        credentials.to_csv(f'{RESULTS_DIR}/passwords.pw', index=False)
    except Exception as e:
        st.error(f"Error adding new user: {e}")


def reset_state():
    """
    Reset the Streamlit session state except for essential keys.
    """
    keys_to_keep = ['logged_in', 'username']
    for key in list(st.session_state.keys()):
        if key not in keys_to_keep:
            del st.session_state[key]


def main():
    """
    The main function handling the overall application workflow, including:
    - User authentication
    - Project selection/creation
    - Docking process setup and submission
    - Viewing and downloading results
    """
    # Ensure static directory exists
    os.makedirs('static', exist_ok=True)

    # Start shared HTTP server (only once)
    start_shared_http_server()

    # Initialize session state variables
    if 'logged_in' not in st.session_state:
        st.session_state.logged_in = False
    if 'username' not in st.session_state:
        st.session_state.username = ''
    if 'module' not in st.session_state:
        st.session_state.module = ''
    if 'progress' not in st.session_state:
        st.session_state.progress = 0
    if 'show_users' not in st.session_state:
        st.session_state.show_users = False

    # Authentication section
    if not st.session_state.logged_in:
        if 'show_register' in st.session_state and st.session_state.show_register:
            st.title("Register New User")
            st.write("Please enter a username and password.")
            new_username = st.text_input("New Username")
            new_password = st.text_input("New Password", type='password')
            confirm_password = st.text_input("Confirm Password", type='password')

            if st.button("Register"):
                if not re.match('^[a-z]+$', new_username):
                    st.error("Username must consist of lowercase letters only.")
                elif new_password != confirm_password:
                    st.error("Passwords do not match.")
                elif username_exists(new_username):
                    st.error("Username already exists.")
                else:
                    hashed_password = hash_password(new_password)
                    add_new_user(new_username, hashed_password)
                    st.success("User registered successfully. You can now log in.")
                    st.session_state.registration_successful = True

            if st.button("Return to Log in", key="return_to_login"):
                st.session_state.clear()
                st.rerun()

        else:
            st.title("Auto Dock Vina web-workflow")
            st.write("Please enter your username and password.")
            username_input = st.text_input("Username")
            password_input = st.text_input("Password", type='password')

            if st.button("Login"):
                if check_credentials(username_input, password_input):
                    user_static_folder = os.path.join('static', username_input)
                    if os.path.exists(user_static_folder):
                        shutil.rmtree(user_static_folder)
                        st.info(f"Removed existing static folder for user '{username_input}'")
                    os.makedirs(user_static_folder, exist_ok=True)
                    st.session_state.logged_in = True
                    st.session_state.username = username_input
                    st.rerun()
                else:
                    st.error("Wrong username or password. Please try again.")

            # Note: Informational text block with instructions
            st.markdown(
                '<div style="max-width: 700px; border: 1px solid; padding: 10px;">'
                'If you want to add a new user, remember that the login can only consist of lowercase letters without additional characters or numbers. The password can contain only basic letters and numbers.'
                '</div>'
                '</br>',
                unsafe_allow_html=True
            )

            if st.button("Add New User"):
                st.session_state.show_register = True
                st.rerun()

            if st.button("SHOW USERS"):
                st.session_state.show_users = not st.session_state.show_users
            if st.session_state.show_users:
                try:
                    credentials = pd.read_csv(f'{RESULTS_DIR}/passwords.pw')
                    users_list = credentials['user'].tolist()
                    st.write("Available Users:")
                    st.write(users_list)
                except Exception as e:
                    st.error(f"Error reading users: {e}")

            st.markdown("""
---
<div style="font-size: 12px; color: gray; text-align: center;">
    Docking system based on <b>AutoDock Vina v.1.2.5</b> and <b>P2RANK v.2.4.2</b><br>
    Author: <b>Arkadiusz Leniak</b> email: arkadiusz.leniak@gmail.com<br>
    <a href="https://github.com/Prospero1988/AutoDock_vina_pipeline" target="_blank">GitHub</a>
</div>
""", unsafe_allow_html=True)

    else:
        # If authenticated, display main menu or selected module
        if st.session_state.module == '':
            st.title(f"Welcome, {st.session_state.username}!")
            st.write("Please select a module to continue:")

            modules = [
                'DOCKING', 'QUEUE', 'SHOW RESULTS',
                'DOWNLOAD RESULTS', 'DELETE RESULTS',
                'PyMOL Installation GUIDE', 'LOG OUT'
            ]
            keys = [
                'docking', 'queue', 'show_results',
                'download', 'delete', 'install_guide', 'logout'
            ]

            for module_name, key in zip(modules, keys):
                if st.button(module_name, key=key):
                    if module_name == 'LOG OUT':
                        user_static_folder = os.path.join('static', st.session_state.username)
                        if os.path.exists(user_static_folder):
                            shutil.rmtree(user_static_folder)
                            st.info(f"Removed static folder for user '{st.session_state.username}'")
                        reset_state()
                        st.session_state.logged_in = False
                        st.rerun()
                    else:
                        st.session_state.module = module_name
                        if module_name == 'DOCKING':
                            st.session_state.progress = 1
                        st.rerun()
                st.write(" ")

        elif st.session_state.module == 'DOCKING':
            docking_module()
        elif st.session_state.module == 'QUEUE':
            queue_module()
        elif st.session_state.module == 'DOWNLOAD RESULTS':
            download_results_module()
        elif st.session_state.module == 'DELETE RESULTS':
            delete_results_module()
        elif st.session_state.module == 'SHOW RESULTS':
            results_module()
        elif st.session_state.module == 'PyMOL Installation GUIDE':
            install_guide_module()


def docking_module():
    """
    Handle the docking workflow steps:
    1. Project setup
    2. PDB codes input and chain selection
    3. Ligand file upload
    4. Docking parameters configuration
    5. Summary and job submission
    """
    st.title("DOCKING Module")

    # Return to main menu
    st.markdown("<hr>", unsafe_allow_html=True)
    if st.button("Return to MENU", key='return_to_menu_docking'):
        reset_state()
        st.rerun()

    progress = st.session_state.progress

    # Initialize variables required for the docking steps
    if 'project_name' not in st.session_state:
        st.session_state.project_name = ''
    if 'project_valid' not in st.session_state:
        st.session_state.project_valid = False
    if 'project_exists' not in st.session_state:
        st.session_state.project_exists = False
    if 'confirmed_project_creation' not in st.session_state:
        st.session_state.confirmed_project_creation = False

    # STEP 1: Project Setup
    if progress == 1:
        st.header("1. Project Setup")
        st.write("Provide a project name for your docking results.")
        project_name_input = st.text_input("Project Name", value=st.session_state.project_name)

        if st.button("Submit Project Name", key="submit_project_name"):
            if project_name_input:
                project_name = validate_project_name(project_name_input)
                if project_name:
                    prefixed_project_name = f"{st.session_state.username}_{project_name}"
                    project_path = os.path.join(RESULTS_DIR, prefixed_project_name)
                    template_path = os.path.join(RESULTS_DIR, 'template')

                    st.session_state.project_name = project_name
                    st.session_state.prefixed_project_name = prefixed_project_name

                    if not os.path.exists(project_path):
                        # Create a new project from the template
                        try:
                            shutil.copytree(template_path, project_path)
                            st.success(f"Project folder '{prefixed_project_name}' has been created.")
                            st.session_state.project_just_created = True
                        except Exception as e:
                            st.error(f"Error copying template: {e}")
                    else:
                        # Project with given name exists
                        st.warning(f"Project '{project_name}' already exists.")
                        st.session_state.project_exists = True
                        st.rerun()
                else:
                    st.error("Invalid project name. Use only alphanumeric characters, numbers, or underscores (_).")
                    st.session_state.project_valid = False
            else:
                st.error("Please enter a project name.")

        if 'project_just_created' in st.session_state and st.session_state.project_just_created and progress == 1:
            if st.button("CONFIRM", key="confirm_project_creation"):
                st.session_state.project_valid = True
                st.session_state.project_exists = False
                st.session_state.project_just_created = False
                st.session_state.progress = 2
                st.rerun()

    if progress == 1 and st.session_state.project_exists:
        st.warning(f"⚠ Project '{st.session_state.project_name}' already exists.")
        decision = st.radio("Project already exists. Choose an action:", ("PROCEED", "CHANGE NAME"), key="decision_action")
        if st.button("Confirm Action", key="confirm_action"):
            if decision == "PROCEED":
                st.session_state.project_valid = True
                st.session_state.progress = 2
                st.session_state.project_exists = False
                st.rerun()
            elif decision == "CHANGE NAME":
                st.session_state.project_valid = False
                st.session_state.project_name = ''
                st.session_state.project_exists = False
                st.session_state.progress = 1
                st.rerun()

    # STEP 2: PDB Codes Input
    if st.session_state.project_valid and st.session_state.progress == 2:
        st.header("2. Enter PDB Codes")
        st.write("Enter PDB codes separated by commas or upload a CSV with PDB codes and chains.")
        pdb_input = st.text_area("PDB Codes", key='pdb_input')
        st.write("Or upload a CSV file with two columns: PDB code and chain ID (no headers).")
        pdb_file = st.file_uploader("Upload CSV", type=['csv'], key='pdb_file')

        project_receptors_path = os.path.join(RESULTS_DIR, st.session_state.prefixed_project_name, 'receptors')
        receptors_csv_path = os.path.join(project_receptors_path, 'receptors.csv')
        pdbs_dir = os.path.join(project_receptors_path, 'pdbs')
        os.makedirs(pdbs_dir, exist_ok=True)

        if 'pdb_codes' not in st.session_state:
            st.session_state.pdb_codes = []
        if 'chains_available' not in st.session_state:
            st.session_state.chains_available = {}
        if 'chains_selected' not in st.session_state:
            st.session_state.chains_selected = {}
        if 'uploaded_pdb_chains' not in st.session_state:
            st.session_state.uploaded_pdb_chains = {}
        if 'download_progress' not in st.session_state:
            st.session_state.download_progress = 0
        if 'total_pdbs' not in st.session_state:
            st.session_state.total_pdbs = 0
        if 'receptors_confirmed' not in st.session_state:
            st.session_state.receptors_confirmed = False

        if st.button("Submit PDB Codes"):
            if pdb_file:
                # Load PDB codes from CSV
                try:
                    df = pd.read_csv(pdb_file, header=None)
                    if df.shape[1] != 2:
                        st.error("CSV must have exactly two columns: PDB code, Chain ID.")
                    else:
                        st.session_state.pdb_codes = df.iloc[:, 0].astype(str).str.strip().str.upper().tolist()
                        st.session_state.uploaded_pdb_chains = dict(zip(
                            df.iloc[:, 0].astype(str).str.strip().str.upper(),
                            df.iloc[:, 1].astype(str).str.strip()
                        ))
                        st.success("CSV file has been successfully loaded.")
                except Exception as e:
                    st.error(f"Error reading CSV file: {e}")
            elif pdb_input:
                # Load PDB codes from manual input
                pdb_input_clean = re.sub(r'\s+', '', pdb_input)
                st.session_state.pdb_codes = [code.strip().upper() for code in pdb_input_clean.split(',') if code.strip()]
                st.session_state.uploaded_pdb_chains = {}
                if st.session_state.pdb_codes:
                    st.success("PDB codes have been loaded from manual entry.")
                else:
                    st.error("Please enter valid PDB codes.")
            else:
                st.error("Please enter PDB codes or upload a CSV file.")

        # Download and parse PDB files to identify chains
        if st.session_state.pdb_codes and not st.session_state.chains_available:
            st.info("Downloading and parsing PDB/mmCIF files...")
            parser = PDBParser(QUIET=True)
            failed_downloads = []
            st.session_state.total_pdbs = len(st.session_state.pdb_codes)
            progress_bar = st.progress(0)

            for idx, pdb_code in enumerate(st.session_state.pdb_codes, start=1):
                pdb_file_path = os.path.join(pdbs_dir, f"{pdb_code}.pdb")
                if not os.path.exists(pdb_file_path):
                    # Attempt downloading PDB first, then mmCIF if PDB fails
                    try:
                        url_pdb = f"https://files.rcsb.org/download/{pdb_code}.pdb"
                        response_pdb = requests.get(url_pdb)
                        if response_pdb.status_code == 200:
                            with open(pdb_file_path, 'w') as f:
                                f.write(response_pdb.text)
                            structure = parser.get_structure(pdb_code, pdb_file_path)
                            chains = sorted([chain.id for chain in structure.get_chains()])
                            st.session_state.chains_available[pdb_code] = chains
                        else:
                            # Try mmCIF format
                            url_cif = f"https://files.rcsb.org/download/{pdb_code}.cif"
                            response_cif = requests.get(url_cif)
                            if response_cif.status_code == 200:
                                cif_file_path = os.path.join(pdbs_dir, f"{pdb_code}.cif")
                                with open(cif_file_path, 'w') as f:
                                    f.write(response_cif.text)

                                # Convert mmCIF to PDB using OBABEL
                                converted_pdb_path = pdb_file_path
                                try:
                                    subprocess.run([
                                        OBABEL_PATH,
                                        "-i", "cif",
                                        cif_file_path,
                                        "-o", "pdb",
                                        "-O", converted_pdb_path
                                    ], check=True)

                                    if os.path.exists(converted_pdb_path):
                                        structure = parser.get_structure(pdb_code, converted_pdb_path)
                                        chains = sorted([chain.id for chain in structure.get_chains()])
                                        st.session_state.chains_available[pdb_code] = chains
                                        os.remove(cif_file_path)
                                    else:
                                        raise FileNotFoundError(f"Conversion failed for {pdb_code}.")
                                except subprocess.CalledProcessError as e:
                                    st.error(f"Error converting mmCIF to PDB for {pdb_code}: {e}")
                                    failed_downloads.append(pdb_code)
                            else:
                                failed_downloads.append(pdb_code)
                    except Exception as e:
                        st.error(f"Error processing PDB code {pdb_code}: {e}")
                        failed_downloads.append(pdb_code)
                else:
                    # If PDB already exists locally
                    try:
                        structure = parser.get_structure(pdb_code, pdb_file_path)
                        chains = sorted([chain.id for chain in structure.get_chains()])
                        st.session_state.chains_available[pdb_code] = chains
                    except Exception as e:
                        st.error(f"Error parsing existing PDB file for {pdb_code}: {e}")
                        failed_downloads.append(pdb_code)

                st.session_state.download_progress = idx / st.session_state.total_pdbs
                progress_bar.progress(st.session_state.download_progress)

            if failed_downloads:
                st.error(f"Failed to download or convert: {', '.join(failed_downloads)}")
            else:
                st.success("All PDB/mmCIF files processed successfully.")
            progress_bar.empty()

        # Chain selection UI
        if st.session_state.chains_available:
            st.header("Select Chain for Each Receptor")
            for pdb_code in st.session_state.pdb_codes:
                if pdb_code in st.session_state.chains_available:
                    available_chains = st.session_state.chains_available[pdb_code]
                    if not available_chains:
                        st.warning(f"No chains found for {pdb_code}.")
                        continue
                    desired_chain = st.session_state.uploaded_pdb_chains.get(pdb_code, available_chains[0])
                    if desired_chain not in available_chains:
                        selected_chain = available_chains[0]
                        st.warning(f"Chain '{desired_chain}' not found. Using '{selected_chain}'.")
                    else:
                        selected_chain = desired_chain

                    selected = st.selectbox(
                        f"Select chain for {pdb_code}",
                        options=available_chains,
                        index=available_chains.index(selected_chain),
                        key=f"chain_select_{pdb_code}"
                    )
                    st.session_state.chains_selected[pdb_code] = selected

            if st.button("Confirm Receptors"):
                if st.session_state.chains_selected:
                    try:
                        receptors_df = pd.DataFrame([
                            [pdb, chain] for pdb, chain in st.session_state.chains_selected.items()
                        ])
                        receptors_df.to_csv(receptors_csv_path, index=False, header=False)
                        st.success(f"Receptors saved in {os.path.basename(receptors_csv_path)}.")
                        st.session_state.receptors_confirmed = True
                    except Exception as e:
                        st.error(f"Error saving receptors.csv: {e}")
                else:
                    st.error("No chains selected.")

            if st.session_state.receptors_confirmed:
                if st.button("Proceed to Next Step"):
                    try:
                        if os.path.exists(pdbs_dir):
                            shutil.rmtree(pdbs_dir)
                            st.success("Temporary PDB directory removed.")
                        st.session_state.progress = 3
                        st.rerun()
                    except Exception as e:
                        st.error(f"Error removing temporary directory: {e}")

    # STEP 3: Upload Ligand
    if st.session_state.project_valid and st.session_state.progress == 3:
        st.header("3. Upload Ligand Files")
        ligand_file = st.file_uploader("Choose a ligand file (.mol2 or .SDF)", type=['mol2', 'SDF', 'sdf'], key='ligand_file')

        ligands_folder = os.path.join(RESULTS_DIR, st.session_state.prefixed_project_name, 'ligands')
        os.makedirs(ligands_folder, exist_ok=True)

        if st.button("Upload Ligand File"):
            if ligand_file:
                ligand_file_name = ligand_file.name
                if ligand_file_name.lower().endswith(('.mol2', '.sdf')):
                    ligand_file_path = os.path.join(ligands_folder, ligand_file_name)
                    with open(ligand_file_path, 'wb') as f:
                        f.write(ligand_file.getbuffer())
                    st.success("Ligand file uploaded.")
                    st.session_state.ligand_file_name = ligand_file_name
                    st.session_state.ligand_uploaded = True
                    st.session_state.progress = 4
                    st.rerun()
                else:
                    st.error("Invalid file format. Only .mol2 and .SDF files accepted.")
            else:
                st.error("Please upload a ligand file.")

    # STEP 4: Docking Parameters
    if st.session_state.project_valid and st.session_state.progress == 4:
        st.header("4. Docking Parameters")
        st.write("Adjust docking parameters if needed, or proceed with defaults.")

        parameters = {
            'pckt': {'value': '1', 'default': '1', 'description': 'Pocket number (default: 1).'},
            'exhaust': {'value': '16', 'default': '16', 'description': 'Exhaustiveness (default: 16).'},
            'energy_range': {'value': '4', 'default': '4', 'description': 'Energy range (default: 4).'},
            'tol_x': {'value': '', 'default': '', 'description': 'Tolerance in X dimension.'},
            'tol_y': {'value': '', 'default': '', 'description': 'Tolerance in Y dimension.'},
            'tol_z': {'value': '', 'default': '', 'description': 'Tolerance in Z dimension.'},
            'offset_x': {'value': '0', 'default': '0', 'description': 'Offset in X (default: 0).'},
            'offset_y': {'value': '0', 'default': '0', 'description': 'Offset in Y (default: 0).'},
            'offset_z': {'value': '0', 'default': '0', 'description': 'Offset in Z (default: 0).'},
        }

        if 'parameters_set' not in st.session_state:
            st.session_state.parameters_set = {}

        for param in parameters:
            col1, col2 = st.columns([1, 3])
            with col1:
                checkbox = st.checkbox(param, value=(param in st.session_state.parameters_set), help=parameters[param]['description'], key=f"{param}_checkbox")
            with col2:
                if checkbox:
                    value = st.text_input(
                        f"Value for {param}",
                        value=st.session_state.parameters_set.get(param, parameters[param]['default']),
                        key=param,
                        help=parameters[param]['description']
                    )
                    parameters[param]['value'] = value
                    st.session_state.parameters_set[param] = value
                else:
                    parameters[param]['value'] = ''
                    if param in st.session_state.parameters_set:
                        del st.session_state.parameters_set[param]

        if st.button("Proceed to Summary"):
            st.session_state.progress = 5
            st.session_state.parameters = parameters
            st.rerun()

    # STEP 5: Summary
    if st.session_state.project_valid and st.session_state.progress == 5:
        st.header("5. Project Summary")
        pdb_codes = []
        receptors_csv_path = os.path.join(RESULTS_DIR, st.session_state.prefixed_project_name, 'receptors', 'receptors.csv')
        if os.path.exists(receptors_csv_path):
            with open(receptors_csv_path, 'r') as f:
                pdb_codes = [line.strip() for line in f if line.strip()]

        parameters = st.session_state.parameters
        st.write(f"**Project Name:** {st.session_state.project_name}")
        st.write(f"**PDB Codes:** {', '.join(pdb_codes) if pdb_codes else 'None'}")
        st.write(f"**Ligand File:** {st.session_state.ligand_file_name if 'ligand_file_name' in st.session_state else 'None'}")
        st.write("**Docking Parameters:**")
        for param in parameters:
            if parameters[param]['value']:
                st.write(f"- {param}: {parameters[param]['value']}")
            else:
                st.write(f"- {param}: (Not set)")

        if st.button("Start Docking"):
            st.session_state.progress = 6
            st.rerun()

    # STEP 6: Start Docking
    if st.session_state.project_valid and st.session_state.progress == 6:
        st.header("6. Start Docking")
        if 'ligand_file_name' not in st.session_state or not st.session_state.ligand_file_name:
            st.error("No ligand file uploaded.")
        else:
            pdb_codes = []
            receptors_csv_path = os.path.join(RESULTS_DIR, st.session_state.prefixed_project_name, 'receptors', 'receptors.csv')
            if os.path.exists(receptors_csv_path):
                with open(receptors_csv_path, 'r') as f:
                    pdb_codes = [line.strip() for line in f if line.strip()]

            if not pdb_codes:
                st.error("No PDB codes provided.")
            else:
                parameters = st.session_state.parameters
                script_content = f"""#!/bin/bash
#SBATCH --job-name={st.session_state.username}_{st.session_state.project_name}
#SBATCH --output=docking_output.log
#SBATCH --error=docking_error.log
#SBATCH --ntasks=1
#SBATCH --time=INFINITE
#SBATCH --partition=main

source ~/miniconda/etc/profile.d/conda.sh
conda activate auto_dock

python3 init_docking.py --pdb_ids receptors.csv --ligands '{st.session_state.ligand_file_name}'"""

                for param in parameters:
                    if parameters[param]['value']:
                        script_content += f" --{param} {parameters[param]['value']}"

                script_content += f"\necho 'Job submitted by {st.session_state.username}'\n"

                script_path = os.path.join(RESULTS_DIR, st.session_state.prefixed_project_name, 'start_docking.sh')
                with open(script_path, 'w') as f:
                    f.write(script_content)
                os.chmod(script_path, 0o755)
                try:
                    subprocess.run(['sbatch', script_path], cwd=os.path.dirname(script_path))
                    st.success("Docking job has been submitted to the queue.")
                    st.session_state.progress = 1
                except Exception as e:
                    st.error(f"Error submitting job: {e}")


def queue_module():
    """
    Display the SLURM queue and allow the user to cancel their own jobs.
    """
    st.title("QUEUE Module")
    st.write("Current SLURM queue:")

    def display_queue():
        try:
            cmd = ['squeue', '-r', '-o', '%i,%j,%T,%M,%V', '--noheader']
            result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
            output = result.stdout.strip().split('\n')
            rows = [line.strip().split(',') for line in output if line.strip()]
            if rows:
                df = pd.DataFrame(rows, columns=['JobID', 'JobName', 'State', 'TimeUsed', 'SubmitTime'])
                st.table(df)
            else:
                st.write("No jobs in the queue.")
        except Exception as e:
            st.error(f"Error fetching SLURM queue: {e}")

    display_queue()

    if st.button("REFRESH"):
        st.rerun()

    st.write("Cancel a job by providing its JOB_ID:")
    job_id_input = st.text_input("Enter JOB_ID")

    if st.button("Cancel Job"):
        cmd_check = ['squeue', '-j', job_id_input, '-o', '%j', '--noheader']
        result_check = subprocess.run(cmd_check, stdout=subprocess.PIPE, text=True)
        job_name = result_check.stdout.strip()

        if job_name:
            if job_name.startswith(f"{st.session_state.username}_"):
                cmd_cancel = ['scancel', job_id_input]
                result_cancel = subprocess.run(cmd_cancel)
                if result_cancel.returncode == 0:
                    st.success(f"Job {job_id_input} has been cancelled.")
                else:
                    st.error(f"Failed to cancel job {job_id_input}.")
            else:
                st.error("This job does not belong to you.")
        else:
            st.error("No job found with the given JOB_ID.")

    if st.button("Return to MENU", key='return_to_menu_queue'):
        reset_state()
        st.rerun()


def download_results_module():
    """
    Allow the user to select their projects and download them as a ZIP archive.
    """
    st.title("DOWNLOAD RESULTS Module")
    dock_folder = RESULTS_DIR
    user_prefix = f"{st.session_state.username}_"
    user_projects = [f for f in os.listdir(dock_folder) if f.startswith(user_prefix)]
    project_names = [f.replace(user_prefix, '') for f in user_projects]

    if not project_names:
        st.info("You have no projects to display.")
        if st.button("Return to MENU", key='return_to_menu_download'):
            reset_state()
            st.rerun()
        return

    st.write("Select projects to download:")
    selected_projects = st.multiselect("Your Projects", project_names, key='selected_projects_download')

    if selected_projects:
        if st.button("Download selected projects"):
            zip_filename = f"docking_results_{st.session_state.username}.zip"
            zip_path = os.path.join('/tmp', zip_filename)

            with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
                for proj in selected_projects:
                    proj_folder = os.path.join(dock_folder, user_prefix + proj)
                    for foldername, subfolders, filenames in os.walk(proj_folder):
                        for filename in filenames:
                            file_path = os.path.join(foldername, filename)
                            arcname = os.path.relpath(file_path, dock_folder)
                            zipf.write(file_path, arcname)

            with open(zip_path, 'rb') as f:
                st.download_button(
                    label="Download ZIP",
                    data=f,
                    file_name=zip_filename,
                    mime='application/zip'
                )

    st.markdown("<hr>", unsafe_allow_html=True)
    if st.button("Return to MENU", key='return_to_menu_download'):
        reset_state()
        st.rerun()


def delete_results_module():
    """
    Allow the user to delete selected projects from their results directory.
    """
    st.title("DELETE RESULTS Module")
    dock_folder = RESULTS_DIR
    user_prefix = f"{st.session_state.username}_"
    user_projects = [f for f in os.listdir(dock_folder) if f.startswith(user_prefix)]
    project_names = [f.replace(user_prefix, '') for f in user_projects]

    if not project_names:
        st.info("You have no projects to display.")
        if st.button("Return to MENU", key='return_to_menu_delete'):
            reset_state()
            st.rerun()
        return

    st.write("Select projects to delete:")
    selected_projects = st.multiselect("Your Projects", project_names, key='selected_projects_delete')

    if selected_projects:
        if 'confirm_delete' not in st.session_state:
            if st.button("Delete selected projects"):
                st.session_state.confirm_delete = True
                st.rerun()
        else:
            if st.button("Confirm Deletion", key='confirm_deletion'):
                for proj in selected_projects:
                    proj_folder = os.path.join(dock_folder, user_prefix + proj)
                    try:
                        shutil.rmtree(proj_folder)
                        st.success(f"Project '{proj}' has been deleted.")
                    except Exception as e:
                        st.error(f"Error deleting '{proj}': {e}")
                del st.session_state.confirm_delete
                st.rerun()
            else:
                st.warning("Press 'Confirm Deletion' to permanently delete selected projects.")
    else:
        st.info("No projects selected.")

    st.markdown("<hr>", unsafe_allow_html=True)
    if st.button("Return to MENU", key='return_to_menu_delete'):
        reset_state()
        st.rerun()


def results_module():
    """
    Display results for a selected project and receptor, allowing the user to:
    - View interactive HTML results
    - Download CSV result files
    """
    st.title("SHOW RESULTS Module")

    if 'host_ip' not in st.session_state:
        st.session_state.host_ip = 'localhost'

    st.write("If connecting locally, leave 'localhost'. Otherwise, specify the host IP/domain.")
    host_input = st.text_input("Host address (IP or domain):", value=st.session_state.host_ip)

    if st.button("CONFIRM HOST"):
        st.session_state.host_ip = host_input
        st.success(f"Host address set to: {st.session_state.host_ip}")

    dock_folder = RESULTS_DIR
    user_prefix = f"{st.session_state.username}_"
    user_projects = [f for f in os.listdir(dock_folder) if f.startswith(user_prefix)]
    project_names = [f.replace(user_prefix, '') for f in user_projects]

    if not project_names:
        st.info("You have no projects to display.")
        if st.button("Return to MENU", key='return_to_menu_results'):
            reset_state()
            st.rerun()
        return

    st.write("Select a project:")
    selected_project = st.selectbox("Your Projects", project_names, key='selected_project')

    if selected_project:
        project_folder = os.path.join(dock_folder, user_prefix + selected_project)
        receptors_csv_path = os.path.join(project_folder, 'receptors', 'receptors.csv')

        if os.path.exists(receptors_csv_path):
            df_receptors = pd.read_csv(receptors_csv_path, header=None, names=['PDB_ID', 'Chain_ID'])
            df_receptors['PDB_ID'] = df_receptors['PDB_ID'].str.strip().str.upper()
            df_receptors['Chain_ID'] = df_receptors['Chain_ID'].str.strip().str.upper()
            df_receptors['Combined'] = df_receptors['PDB_ID'] + "_" + df_receptors['Chain_ID']
            receptors = df_receptors['Combined'].tolist()
        else:
            st.error(f"Receptors file not found: {receptors_csv_path}")
            return

        if receptors:
            st.write("Select a receptor:")
            selected_receptor = st.selectbox("Receptors", receptors, key='selected_receptor')

            if selected_receptor:
                receptor_folder = os.path.join(project_folder, selected_receptor)
                if os.path.exists(receptor_folder):
                    html_files = [f for f in os.listdir(receptor_folder) if f.endswith('.html')]
                    if html_files:
                        html_file_name = html_files[0]

                        if st.button("SHOW INTERACTIVE RESULTS"):
                            static_receptor_folder = os.path.join('static', st.session_state.username, selected_project, selected_receptor)
                            os.makedirs(os.path.dirname(static_receptor_folder), exist_ok=True)
                            if os.path.exists(static_receptor_folder):
                                shutil.rmtree(static_receptor_folder)

                            try:
                                shutil.copytree(receptor_folder, static_receptor_folder)
                            except Exception as e:
                                st.error(f"Error copying files: {e}")
                                return

                            html_url = f"http://{st.session_state.host_ip}:{HTTP_SERVER_PORT}/{st.session_state.username}/{selected_project}/{selected_receptor}/{html_file_name}"
                            st.markdown(f'<a href="{html_url}" target="_blank">Open Interactive Results in New Tab</a>', unsafe_allow_html=True)
                    else:
                        st.error("No HTML file found for this receptor.")

                    st.markdown("<hr>", unsafe_allow_html=True)
                    csv_files = [f for f in os.listdir(receptor_folder) if f.endswith('.csv')]
                    if csv_files:
                        csv_file_path = os.path.join(receptor_folder, csv_files[0])
                        with open(csv_file_path, 'rb') as f:
                            csv_data = f.read()
                        st.download_button(
                            label="DOWNLOAD RESULTS IN CSV",
                            data=csv_data,
                            file_name=csv_files[0],
                            mime='text/csv'
                        )
                    else:
                        st.error("No CSV file found for this receptor.")
                else:
                    st.error(f"Receptor folder '{selected_receptor}' does not exist.")
        else:
            st.error("No receptors found in receptors.csv.")

    if st.button("Return to MENU", key='return_to_menu_results'):
        reset_state()
        st.rerun()


def install_guide_module():
    """
    Display the installation guide for Open-PyMOL on Windows.
    """
    st.title("Installation Guide for Open-PyMOL on Windows")

    html_content = """
    <div class="frame">
        <h2>1. Download and Install Miniconda</h2>
        <ol>
            <li>Navigate to the <a href="https://docs.conda.io/en/latest/miniconda.html" target="_blank">Miniconda download page</a>.</li>
            <li>Download the appropriate Miniconda installer for Windows.</li>
            <li>Run the installer and follow the instructions.</li>
            <li>Verify installation:
                <pre><code>conda --version</code></pre>
            </li>
        </ol>

        <h2>2. Initialize Conda, Create and Activate the Open_PyMOL Environment</h2>
        <ol>
            <li>Open the Command Prompt.</li>
            <li>
                <pre><code>conda init</code></pre>
            </li>
            <li>Restart the Command Prompt.</li>
            <li>Create a new environment:
                <pre><code>conda create -n open_pymol python=3.10</code></pre>
            </li>
            <li>Activate the environment:
                <pre><code>conda activate open_pymol</code></pre>
            </li>
        </ol>

        <h2>3. Download and Install Open-PyMOL</h2>
        <ol>
            <li>
                <pre><code>conda install -y -c conda-forge pymol-open-source</code></pre>
            </li>
        </ol>

        <h2>Running Open-PyMOL</h2>
        <ol>
            <li>Open Command Prompt.</li>
            <li><pre><code>conda activate open_pymol</code></pre></li>
            <li><pre><code>pymol</code></pre></li>
        </ol>
    </div>
    """

    frame_style = """
    <style>
        .frame {
            max-width: 1600px;
            width: 100%;
            background-color: white;
            color: black;
            padding: 20px;
            box-sizing: border-box;
            overflow-x: auto;
            margin: 0 auto;
        }
        .frame h2 {
            color: black;
        }
        .frame a {
            color: #1E90FF;
            text-decoration: none;
        }
        .frame a:hover {
            text-decoration: underline;
        }
        .frame ol {
            margin-left: 20px;
        }
        .frame pre {
            background-color: #f4f4f4;
            padding: 10px;
            border-radius: 4px;
            overflow-x: auto;
        }
        .frame code {
            font-family: Consolas, "Courier New", monospace;
            color: #c7254e;
            background-color: #f9f2f4;
            padding: 2px 4px;
            border-radius: 4px;
        }
    </style>
    """

    components.html(frame_style + html_content, height=1150, scrolling=True)

    st.markdown("<hr>", unsafe_allow_html=True)
    if st.button("RETURN TO MENU", key='return_to_menu_install_guide'):
        reset_state()
        st.rerun()


if __name__ == "__main__":
    main()