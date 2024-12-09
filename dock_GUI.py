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

# Set page configuration
st.set_page_config(page_title="Docking Program", layout="centered")

# Global settings for HTTP server
HTTP_SERVER_PORT = 8001  # You can change this port if needed
HTTPD = None  # Global HTTP Server instance

def start_http_server(directory, port=HTTP_SERVER_PORT):
    global HTTPD
    handler = partial(SimpleHTTPRequestHandler, directory=directory)
    HTTPD = HTTPServer(('0.0.0.0', port), handler)
    HTTPD.serve_forever()

def stop_http_server():
    global HTTPD
    if HTTPD is not None:
        HTTPD.shutdown()
        HTTPD = None

# Ensure 'static' directory exists and clear 'static' when app starts
if os.path.exists('static'):
    shutil.rmtree('static')
os.makedirs('static')

def validate_project_name(name):
    # Replace spaces with underscores
    name = name.replace(' ', '_')
    # Check if name contains only allowed characters
    if re.match(r'^[A-Za-z0-9_]+$', name):
        return name
    else:
        return None

def hash_password(password):
    # Generate salt and hash the password
    salt = bcrypt.gensalt()
    hashed = bcrypt.hashpw(password.encode('utf-8'), salt)
    return hashed.decode('utf-8')  # Store as string

def check_credentials(username, password):
    try:
        credentials = pd.read_csv('passwords.pw')
        user_row = credentials[credentials['user'] == username]
        if not user_row.empty:
            stored_hash = user_row.iloc[0]['password']
            # Check password
            if bcrypt.checkpw(password.encode('utf-8'), stored_hash.encode('utf-8')):
                return True
    except Exception as e:
        print(e)
        return False
    return False

def username_exists(username):
    try:
        credentials = pd.read_csv('passwords.pw')
        return username in credentials['user'].values
    except Exception:
        return False

def add_new_user(username, hashed_password):
    try:
        # Check if passwords.pw exists
        if os.path.exists('passwords.pw'):
            credentials = pd.read_csv('passwords.pw')
            # Append new user
            new_user = pd.DataFrame({'user': [username], 'password': [hashed_password]})
            credentials = pd.concat([credentials, new_user], ignore_index=True)
        else:
            # Create new DataFrame
            credentials = pd.DataFrame({'user': [username], 'password': [hashed_password]})
        # Save to passwords.pw
        credentials.to_csv('passwords.pw', index=False)
    except Exception as e:
        st.error(f"Error adding new user: {e}")

def main():
    # Initialize session state variables
    if 'logged_in' not in st.session_state:
        st.session_state.logged_in = False
    if 'username' not in st.session_state:
        st.session_state.username = ''
    if 'module' not in st.session_state:
        st.session_state.module = ''
    if 'progress' not in st.session_state:
        st.session_state.progress = 0
    if 'http_server_thread' not in st.session_state:
        st.session_state['http_server_thread'] = None

    # Login screen
    if not st.session_state.logged_in:
        if 'show_register' in st.session_state and st.session_state.show_register:
            st.title("Register New User")
            st.write("Please enter a username and password.")

            new_username = st.text_input("New Username")
            new_password = st.text_input("New Password", type='password')
            confirm_password = st.text_input("Confirm Password", type='password')

            if st.button("Register"):
                # Validate username
                if not re.match('^[a-z]+$', new_username):
                    st.error("Username must consist of lowercase letters only.")
                elif new_password != confirm_password:
                    st.error("Passwords do not match.")
                elif username_exists(new_username):
                    st.error("Username already exists.")
                else:
                    # Hash password and store new user
                    hashed_password = hash_password(new_password)
                    add_new_user(new_username, hashed_password)
                    st.success("User registered successfully. You can now log in.")
                    st.session_state.registration_successful = True

            st.write(" ")  # Add more space        

            # Jeden przycisk "Return to Log in"
            if st.button("Return to Log in", key="return_to_login"):
                reset_state()
                st.rerun()

        else:
            st.title("Docking Program Login")
            st.write("Please enter your username and password.")

            username_input = st.text_input("Username")
            password_input = st.text_input("Password", type='password')

            if st.button("Login"):
                if check_credentials(username_input, password_input):
                    # Clear 'static' when any user logs in
                    if os.path.exists('static'):
                        shutil.rmtree('static')
                    os.makedirs('static')
                    st.session_state.logged_in = True
                    st.session_state.username = username_input
                    st.rerun()
                else:
                    st.error("Wrong username or password. Please try again.")

            st.markdown(
                '<div style="max-width: 300px; border: 1px solid; padding: 10px;">'
                'If you want to add a new user, remember that the login can only consist of lowercase letters without additional characters or numbers. The password can contain only basic letters and numbers.'
                '</div>',
                unsafe_allow_html=True
            )
            st.write(" ")  # Add some space

            # Add 'Add New User' and 'SHOW USERS' buttons below the text
            col1, col2 = st.columns(2)
            with col1:
                if st.button("Add New User"):
                    st.session_state.show_register = True
                    st.rerun()
            with col2:
                if st.button("SHOW USERS"):
                    try:
                        credentials = pd.read_csv('passwords.pw')
                        users_list = credentials['user'].tolist()
                        st.write("Available Users:")
                        st.write(users_list)
                    except Exception as e:
                        st.error(f"Error reading users: {e}")

    else:
        # Main Menu
        if st.session_state.module == '':
            st.title(f"Welcome, {st.session_state.username}!")
            st.write("Please select a module to continue:")

            modules = ['DOCKING', 'QUEUE', 'SHOW RESULTS', 'DOWNLOAD RESULTS', 'DELETE RESULTS', 'PyMOL Installation GUIDE', 'LOG OUT']
            keys = ['docking', 'queue', 'show_results', 'download', 'delete', 'install_guide', 'logout']


            st.write(" ")  # Add some space
            for module_name, key in zip(modules, keys):
                if st.button(module_name, key=key):
                    if module_name == 'LOG OUT':
                        # Clean up static directory
                        if os.path.exists('static'):
                            shutil.rmtree('static')
                        os.makedirs('static')
                        # Stop the HTTP server if running
                        stop_http_server()
                        reset_state()
                        st.session_state.logged_in = False
                        st.rerun()
                    else:
                        st.session_state.module = module_name
                        if module_name == 'DOCKING':
                            st.session_state.progress = 1
                        st.rerun()
                st.write(" ")  # Add some space

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
    st.title("DOCKING Module")

    # Add 'Return to MENU' button at the bottom
    st.markdown("<hr>", unsafe_allow_html=True)
    if st.button("Return to MENU", key='return_to_menu_docking'):
        reset_state()
        st.rerun()

    # Proceed through the steps based on progress
    progress = st.session_state.progress

    # Initialize parameters outside the conditional blocks
    parameters = {}

    # Step 1: Greeting and Project Name Input
    if progress == 1:
        st.header("1. Project Setup")
        st.write("Please provide a project name where your docking will be saved.")

        # Initialize variables
        if 'project_name' not in st.session_state:
            st.session_state.project_name = ''
        if 'project_valid' not in st.session_state:
            st.session_state.project_valid = False
        if 'project_exists' not in st.session_state:
            st.session_state.project_exists = False

        # Project Name Input
        project_name_input = st.text_input("Project Name", st.session_state.project_name)

        if st.button("Submit Project Name"):
            if project_name_input:
                project_name = validate_project_name(project_name_input)
                if project_name:
                    st.session_state.project_name = project_name
                    # Add user prefix to project folder name
                    prefixed_project_name = f"{st.session_state.username}_{project_name}"
                    st.session_state.prefixed_project_name = prefixed_project_name
                    project_path = os.path.join('/home/docking_machine/dock', prefixed_project_name)
                    template_path = '/home/docking_machine/dock/template'

                    if not os.path.exists(project_path):
                        # Create the folder
                        try:
                            shutil.copytree(template_path, project_path)
                            st.success(f"Project folder '{prefixed_project_name}' has been created.")
                            st.session_state.project_valid = True
                            st.session_state.progress = 2
                            st.rerun()
                        except Exception as e:
                            st.error(f"Error copying template: {e}")
                    else:
                        # Project already exists
                        st.session_state.project_exists = True
                        st.warning(f"Project '{project_name}' already exists.")
                        use_existing = st.radio(
                            "Do you want to continue with the existing project or enter a new project name?",
                            ('Continue with existing project', 'Enter a new project name')
                        )
                        if use_existing == 'Continue with existing project':
                            st.session_state.project_valid = True
                            st.session_state.progress = 2
                            st.rerun()
                        else:
                            st.session_state.project_valid = False
                            st.session_state.project_name = ''
                else:
                    st.error("Invalid project name. Use only alphanumeric characters, numbers, or underscores (_).")
                    st.session_state.project_valid = False

    # Step 2: Entering PDB codes
    if st.session_state.project_valid and progress == 2:
        st.header("2. Enter PDB Codes")
        st.write("Enter PDB codes, separated by commas (e.g., 1A2B, 2X3C, 3Y4D).")
        pdb_input = st.text_area("PDB Codes", key='pdb_input')
        st.write("Or upload a CSV file with PDB codes (one code per line).")
        pdb_file = st.file_uploader("Upload CSV", type=['csv'], key='pdb_file')

        pdb_codes = []
        project_receptors_path = os.path.join('/home/docking_machine/dock', st.session_state.prefixed_project_name, 'receptors')
        receptors_csv_path = os.path.join(project_receptors_path, 'receptors.csv')
        os.makedirs(project_receptors_path, exist_ok=True)

        if st.button("Submit PDB Codes"):
            if pdb_file:
                # User uploaded a file
                try:
                    df = pd.read_csv(pdb_file, header=None)
                    if df.shape[1] == 1:
                        pdb_codes = df[0].astype(str).tolist()
                        # Save to receptors.csv
                        df.to_csv(receptors_csv_path, index=False, header=False)
                        st.success("PDB codes have been saved.")
                        st.session_state.pdb_codes_saved = True
                        st.session_state.progress = 3
                        st.rerun()
                    else:
                        st.error("Invalid CSV format. Please ensure the file contains one column with PDB codes.")
                except Exception as e:
                    st.error(f"Error reading CSV file: {e}")
            elif pdb_input:
                # User entered PDB codes
                pdb_input_clean = re.sub(r'\s+', '', pdb_input)
                pdb_codes = [code.strip() for code in pdb_input_clean.split(',') if code.strip()]
                if pdb_codes:
                    # Save to receptors.csv
                    with open(receptors_csv_path, 'w') as f:
                        for code in pdb_codes:
                            f.write(f"{code}\n")
                    st.success("PDB codes have been saved.")
                    st.session_state.pdb_codes_saved = True
                    st.session_state.progress = 3
                    st.rerun()
                else:
                    st.error("Please enter valid PDB codes.")
            else:
                st.error("Please enter PDB codes or upload a CSV file.")

    # Step 3: Uploading ligand files
    if st.session_state.project_valid and progress == 3:
        st.header("3. Upload Ligand Files")
        st.write("Upload a ligand file (.mol2 or .SDF).")
        ligand_file = st.file_uploader("Choose a ligand file", type=['mol2', 'SDF', 'sdf'], key='ligand_file')

        ligands_folder = os.path.join('/home/docking_machine/dock', st.session_state.prefixed_project_name, 'ligands')
        os.makedirs(ligands_folder, exist_ok=True)

        if st.button("Upload Ligand File"):
            if ligand_file:
                ligand_file_name = ligand_file.name
                if ligand_file_name.lower().endswith(('.mol2', '.sdf')):
                    ligand_file_path = os.path.join(ligands_folder, ligand_file_name)
                    # Save the uploaded file
                    with open(ligand_file_path, 'wb') as f:
                        f.write(ligand_file.getbuffer())
                    st.success("Ligand file has been uploaded.")
                    st.session_state.ligand_file_name = ligand_file_name
                    st.session_state.ligand_uploaded = True
                    st.session_state.progress = 4
                    st.rerun()
                else:
                    st.error("Invalid file format. Only .mol2 and .SDF files are accepted.")
            else:
                st.error("Please upload a ligand file.")

    # Step 4: Configuration of docking parameters
    if st.session_state.project_valid and progress == 4:
        st.header("4. Docking Parameters")
        st.write("Select and configure docking parameters.")

        # Initialize parameters with default values
        parameters = {
            'tol_x': {'value': '', 'default': '', 'description': 'Tolerance in Ångströms to expand the docking pocket dimension in X beyond those defined by P2Rank (optional).'},
            'tol_y': {'value': '', 'default': '', 'description': 'Tolerance in Ångströms to expand the docking pocket dimension in Y beyond those defined by P2Rank (optional).'},
            'tol_z': {'value': '', 'default': '', 'description': 'Tolerance in Ångströms to expand the docking pocket dimension in Z beyond those defined by P2Rank (optional).'},
            'pckt': {'value': '1', 'default': '1', 'description': 'Pocket number to use from P2Rank predictions (default: 1).'},
            'exhaust': {'value': '16', 'default': '16', 'description': 'Specifies how thorough the search should be for the best binding poses. Higher values increase precision but require more computation time (default: 16).'},
            'energy_range': {'value': '4', 'default': '4', 'description': 'Determines the range of energy scores (in kcal/mol) for poses to be considered (default: 4).'},
            'offset_x': {'value': '0', 'default': '0', 'description': 'Offset in Ångströms to shift the center of the docking grid box along the X-axis (optional, default: 0).'},
            'offset_y': {'value': '0', 'default': '0', 'description': 'Offset in Ångströms to shift the center of the docking grid box along the Y-axis (optional, default: 0).'},
            'offset_z': {'value': '0', 'default': '0', 'description': 'Offset in Ångströms to shift the center of the docking grid box along the Z-axis (optional, default: 0).'},
        }

        if 'parameters_set' not in st.session_state:
            st.session_state.parameters_set = {}

        # For each parameter, display checkbox and input field if checked
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

    # Step 5: Project Summary
    if st.session_state.project_valid and progress == 5:
        st.header("5. Project Summary")
        pdb_codes = []
        receptors_csv_path = os.path.join('/home/docking_machine/dock', st.session_state.prefixed_project_name, 'receptors', 'receptors.csv')
        if os.path.exists(receptors_csv_path):
            with open(receptors_csv_path, 'r') as f:
                pdb_codes = [line.strip() for line in f if line.strip()]

        parameters = st.session_state.parameters

        st.write(f"**Project Name:** {st.session_state.project_name}")
        st.write(f"**PDB Codes:** {', '.join(pdb_codes) if pdb_codes else 'None'}")
        st.write(f"**Ligand File:** {st.session_state.ligand_file_name if st.session_state.ligand_file_name else 'None'}")
        st.write("**Docking Parameters:**")
        for param in parameters:
            if parameters[param]['value']:
                st.write(f"- {param}: {parameters[param]['value']}")
            else:
                st.write(f"- {param}: (Not set)")

        if st.button("Start Docking"):
            st.session_state.progress = 6
            st.rerun()

    # Step 6: Generating and Running Bash Script
    if st.session_state.project_valid and progress == 6:
        st.header("6. Start Docking")
        # Validate that necessary inputs are provided
        if 'ligand_file_name' not in st.session_state or not st.session_state.ligand_file_name:
            st.error("No ligand file uploaded.")
        else:
            pdb_codes = []
            receptors_csv_path = os.path.join('/home/docking_machine/dock', st.session_state.prefixed_project_name, 'receptors', 'receptors.csv')
            if os.path.exists(receptors_csv_path):
                with open(receptors_csv_path, 'r') as f:
                    pdb_codes = [line.strip() for line in f if line.strip()]

            if not pdb_codes:
                st.error("No PDB codes provided.")
            else:
                parameters = st.session_state.parameters
                # Generate start_docking.sh
                # Use the project name with user prefix for job name
                script_content = f"""#!/bin/bash
#SBATCH --job-name={st.session_state.username}_{st.session_state.project_name}
#SBATCH --output=docking_output.log
#SBATCH --error=docking_error.log
#SBATCH --ntasks=1
#SBATCH --time=INFINITE
#SBATCH --partition=main

source ~/miniconda/etc/profile.d/conda.sh
conda activate auto_dock

python3 init_docking.py --pdb_ids receptors.csv --ligands {st.session_state.ligand_file_name}"""

                # Add parameters to the script
                for param in parameters:
                    if parameters[param]['value']:
                        script_content += f" --{param} {parameters[param]['value']}"

                # Optionally, you can add the username to the script environment
                script_content += f"\n\necho 'Job submitted by {st.session_state.username}'"

                # Save the script to start_docking.sh in project folder
                script_path = os.path.join('/home/docking_machine/dock', st.session_state.prefixed_project_name, 'start_docking.sh')
                with open(script_path, 'w') as f:
                    f.write(script_content)
                # Make the script executable
                os.chmod(script_path, 0o755)
                # Run the script using sbatch
                try:
                    subprocess.run(['sbatch', script_path], cwd=os.path.dirname(script_path))
                    st.success("Docking job has been submitted to the queue. Check SLURM logs for progress.")
                    # Reset progress to allow for a new job submission
                    st.session_state.progress = 1
                except Exception as e:
                    st.error(f"Error submitting job: {e}")

def queue_module():
    st.title("QUEUE Module")
    st.write("Current SLURM queue:")

    def display_queue():
        try:
            cmd = ['squeue', '-r', '-o', '%i,%u,%j,%T,%M,%S', '--noheader']
            result = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
            output = result.stdout.strip().split('\n')
            rows = [line.strip().split(',') for line in output if line.strip()]
            if rows:
                df = pd.DataFrame(rows, columns=['JobID', 'User', 'JobName', 'State', 'TimeUsed', 'StartTime'])
                st.table(df)
            else:
                st.write("No jobs in the queue.")
        except Exception as e:
            st.error(f"Error fetching SLURM queue: {e}")

    # Display the queue initially
    display_queue()

    # Add 'REFRESH' button to refresh the queue
    if st.button("REFRESH"):
        st.rerun()

    # Add 'Return to MENU' button at the bottom
    if st.button("Return to MENU", key='return_to_menu_queue'):
        reset_state()
        st.rerun()

def download_results_module():
    st.title("DOWNLOAD RESULTS Module")

    dock_folder = '/home/docking_machine/dock'
    user_prefix = f"{st.session_state.username}_"
    user_projects = [f for f in os.listdir(dock_folder) if f.startswith(user_prefix)]
    project_names = [f.replace(user_prefix, '') for f in user_projects]

    if not project_names:
        st.info("You have no projects to display.")
        if st.button("Return to MENU", key='return_to_menu_download'):
            reset_state()
            st.rerun()
        return

    st.write("Select projects:")
    selected_projects = st.multiselect("Your Projects", project_names, key='selected_projects_download')

    if selected_projects:
        if st.button("Download selected projects"):
            zip_filename = f"docking_results_{st.session_state.username}.zip"
            zip_path = os.path.join('/tmp', zip_filename)
            dock_folder = '/home/docking_machine/dock'
            user_prefix = f"{st.session_state.username}_"

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
    else:
        st.info("No projects selected.")

    st.markdown("<hr>", unsafe_allow_html=True)
    if st.button("Return to MENU", key='return_to_menu_download'):
        reset_state()
        st.rerun()

def delete_results_module():
    st.title("DELETE RESULTS Module")

    dock_folder = '/home/docking_machine/dock'
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
                dock_folder = '/home/docking_machine/dock'
                user_prefix = f"{st.session_state.username}_"
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
    st.title("SHOW RESULTS Module")

    if 'http_server_thread' not in st.session_state:
        st.session_state['http_server_thread'] = None

    dock_folder = '/home/docking_machine/dock'
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
            with open(receptors_csv_path, 'r') as f:
                receptors = [line.strip() for line in f if line.strip()]
        else:
            receptors = []
            st.error("No receptors.csv file found in the project.")

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
                            static_receptor_folder = os.path.join('static', selected_receptor)
                            if os.path.exists(static_receptor_folder):
                                shutil.rmtree(static_receptor_folder)
                            try:
                                shutil.copytree(receptor_folder, static_receptor_folder)
                            except Exception as e:
                                st.error(f"Error copying files: {e}")
                                return

                            # Start HTTP server if not already running
                            if st.session_state['http_server_thread'] is None or not st.session_state['http_server_thread'].is_alive():
                                st.session_state['http_server_thread'] = threading.Thread(target=start_http_server, args=('static',), daemon=True)
                                st.session_state['http_server_thread'].start()

                            html_url = f"http://172.22.31.82:{HTTP_SERVER_PORT}/{selected_receptor}/{html_file_name}"
                            st.markdown(f'<a href="{html_url}" target="_blank">Open Interactive Results in New Tab</a>', unsafe_allow_html=True)
                    else:
                        st.error("No HTML file found in the receptor's folder.")

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
                        st.error("No CSV file found in the receptor's folder.")
                else:
                    st.error(f"Receptor folder '{selected_receptor}' does not exist.")
        else:
            st.error("No receptors found in receptors.csv.")

    if st.button("Return to MENU", key='return_to_menu_results'):
        if os.path.exists('static'):
            shutil.rmtree('static')
        os.makedirs('static')
        reset_state()
        st.rerun()

def install_guide_module():
    st.title("Installation Guide for Open-PyMOL on Windows")

    # Twój kod HTML jako string (bez tagów <html>, <head>, <body>)
    html_content = """
    <div class="frame">
        
        <h2>1. Download and Install Miniconda</h2>
        <ol>
            <li>
                Open your web browser and navigate to the <a href="https://docs.conda.io/en/latest/miniconda.html" target="_blank">Miniconda download page</a>.
            </li>
            <li>
                Download the appropriate Miniconda installer for Windows (64-bit recommended).
            </li>
            <li>
                Once downloaded, run the installer and follow the on-screen instructions. It's recommended to install for "Just Me" and check the option to add Miniconda to your PATH environment variable during installation.
            </li>
            <li>
                After installation, open the Command Prompt to verify the installation by typing:
                <pre><code>conda --version</code></pre>
            </li>
        </ol>

        <h2>2. Initialize Conda, Create and Activate the Open_PyMOL Environment</h2>
        <ol>
            <li>
                Open the Command Prompt.
            </li>
            <li>
                Initialize Conda by running:
                <pre><code>conda init</code></pre>
            </li>
            <li>
                Restart the Command Prompt to apply the changes.
            </li>
            <li>
                Create a new Conda environment named <strong>open_pymol</strong> by executing:
                <pre><code>conda create -n open_pymol python=3.10</code></pre>
            </li>
            <li>
                Activate the newly created environment:
                <pre><code>conda activate open_pymol</code></pre>
            </li>
        </ol>

        <h2>3. Download and Install Open-PyMOL via Conda</h2>
        <ol>
            <li>
                With the <strong>open_pymol</strong> environment activated, install Open-PyMOL by running:
                <pre><code>conda install -y -c conda-forge pymol-open-source</code></pre>
            </li>
            <li>
                Wait for the installation to complete. Conda will handle all necessary dependencies.
            </li>
        </ol>

        <h2>Running Open-PyMOL</h2>
        <ol>
            <li>
                Open the Command Prompt.
            </li>
            <li>
                Activate the <strong>open_pymol</strong> environment:
                <pre><code>conda activate open_pymol</code></pre>
            </li>
            <li>
                Launch Open-PyMOL by typing:
                <pre><code>pymol</code></pre>
            </li>
            <li>
                Open-PyMOL should now start, and you can begin using it for your molecular visualization needs.
            </li>
        </ol>
    </div>
    """

    # Definiowanie stylu ramki
    frame_style = """
    <style>
        .frame {
            max-width: 1600px; /* Zwiększona szerokość o 200px */
            width: 100%; /* Ustawienie szerokości na 100% kontenera */
            background-color: white;
            color: black; /* Ustawienie koloru tekstu na czarny */
            padding: 20px;
            box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);
            box-sizing: border-box;
            overflow-x: auto;
            margin: 0 auto; /* Wycentrowanie ramki */
        }
        .frame h1, .frame h2 {
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

    # Wyświetlanie stylu i HTML w komponentach.html
    components.html(frame_style + html_content, height=1150, scrolling=True)

    st.markdown("<hr>", unsafe_allow_html=True)
    if st.button("RETURN TO MENU", key='return_to_menu_install_guide'):
        reset_state()
        st.rerun()

def reset_state():
    st.session_state.module = ''
    st.session_state.progress = 0
    keys_to_reset = [
        'project_name', 'project_valid', 'project_exists', 'prefixed_project_name',
        'pdb_codes_saved', 'pdb_input', 'pdb_file', 'ligand_file_name',
        'ligand_uploaded', 'parameters_set', 'parameters', 'selected_projects_download',
        'selected_projects_delete', 'confirm_delete', 'show_register', 'registration_successful',
        'selected_project', 'selected_receptor',
        'selected_projects_install_guide'
    ]
    for key in keys_to_reset:
        if key in st.session_state:
            del st.session_state[key]

if __name__ == "__main__":
    main()