import json
import multiprocessing

import streamlit as st
from config import JOB_TRACKER_FILE, REMOTE_HOST, REMOTE_USER
from modules.input_validation import (
    is_unique_id,
    is_valid_email,
    is_valid_smiles,
    is_valid_string,
)
from modules.job_status import JobStatus
from modules.job_submission import JobSubmission
from modules.remote_server import RemoteServer
from modules.utils import (
    load_submitted_jobs,
    remove_empty_fields,
    save_submitted_job,
    zip_job_folder,
)
from rdkit import Chem


def render_tabs(email):
    """
    Renders the tabs in the Streamlit app.
    """
    tab1, tab2, tab3 = st.tabs(
        ["Submit Job", "View Submitted Jobs", "Retrieve Job Output"]
    )

    with tab1:
        render_submit_job_tab(email)

    with tab2:
        render_view_jobs_tab(email)

    with tab3:
        render_retrieve_output_tab(email)


def render_submit_job_tab(email):
    """
    Handles the UI for submitting a job.
    """
    st.header("Submit Job")
    job_name = st.text_input("Job Name")

    if not job_name:
        st.error("Job Name cannot be empty. Please provide a valid job name.")
        return
    valid, error_message = is_valid_string(job_name)
    if not valid:
        st.error(
            f"Invalid Job Name. Only alphabets, numbers, dashes, and underscores are allowed."
        )
        return

    model_seeds = st.text_area("Model Seeds (comma-separated)", value="1")

    unique_ids = set()

    # Collect Protein Sequences
    st.subheader("Define Protein Sequences")
    num_proteins = st.number_input(
        "Number of Protein Entities", min_value=0, value=0, step=1
    )
    protein_sequences = []

    for i in range(num_proteins):
        with st.expander(f"Protein {i + 1} Details"):
            protein_id = st.text_input(
                f"Protein ID (one letter)", key=f"protein_id_{i}"
            )
            sequence = st.text_area(
                f"Protein Sequence", key=f"protein_seq_{i}"
            )

            valid, error_message = is_unique_id(protein_id, unique_ids)
            if not valid:
                st.error(error_message)
            else:
                protein_sequences.append(
                    {
                        "protein": {
                            "id": protein_id.upper(),
                            "sequence": sequence,
                        }
                    }
                )

    # Collect Ligands
    st.subheader("Define Ligands")
    num_ligands = st.number_input(
        "Number of Ligands", min_value=0, value=0, step=1
    )
    ligands = []

    for i in range(num_ligands):
        with st.expander(f"Ligand {i + 1} Details"):
            ligand_id = st.text_input(
                "Ligand ID (one letter)", key=f"ligand_id_{i}"
            )
            smiles = st.text_input(
                "SMILES (optional)", key=f"ligand_smiles_{i}"
            ).strip()

            if not is_valid_smiles(smiles):
                st.error("Invalid SMILES string. Please check your input.")
            else:
                valid, error_message = is_unique_id(ligand_id, unique_ids)
                if not valid:
                    st.error(error_message)
                else:
                    ligands.append(
                        {"ligand": {"id": ligand_id.upper(), "smiles": smiles}}
                    )

    json_content = {
        "name": job_name,
        "modelSeeds": [
            int(seed.strip())
            for seed in model_seeds.split(",")
            if seed.strip()
        ],
        "sequences": protein_sequences + ligands,
        "dialect": "alphafold3",
        "version": 1,
    }

    json_content = remove_empty_fields(json_content)

    st.subheader("Generated JSON")
    st.json(json_content)

    st.download_button(
        label="Download JSON File",
        data=json.dumps(json_content, indent=4),
        file_name=f"{job_name}.json",
        mime="application/json",
    )

    if st.button("Submit Job", key="submit_job"):
        temporary_job_data = {
            "job_id": "Pending",
            "job_name": job_name,
            "status": "Pending",
        }
        save_submitted_job(JOB_TRACKER_FILE, email, temporary_job_data)

        server = RemoteServer(REMOTE_HOST, REMOTE_USER)
        job_submitter = JobSubmission(server)

        process = multiprocessing.Process(
            target=job_submitter.submit, args=(json_content, email, job_name)
        )
        process.start()

        st.success(
            "Your job is being submitted in the background. You can close the app."
        )


def render_view_jobs_tab(email):
    """
    Displays submitted jobs and their statuses.
    """
    st.header("View Submitted Jobs")

    jobs = load_submitted_jobs(JOB_TRACKER_FILE)
    user_jobs = jobs.get(email, [])

    if user_jobs:
        server = RemoteServer(REMOTE_HOST, REMOTE_USER)
        job_status_checker = JobStatus(server)

        st.info("Checking job statuses... Please do not close the app.")
        with st.spinner("Retrieving job statuses, please wait..."):
            for job in user_jobs:
                with st.expander(
                    f"Job: {job['job_name']} (Status: {job['status']})"
                ):
                    st.write(f"**Job Name**: {job['job_name']}")
                    st.write(f"**Job ID**: {job.get('job_id', 'N/A')}")

                    if job["job_id"] == "Failed":
                        st.warning("This job submission failed.")
                        continue

                    current_status, error = job_status_checker.check_status(
                        job["job_id"]
                    )
                    if error:
                        st.error(f"Error checking job status: {error}")
                    else:
                        st.write(f"**Current Status**: {current_status}")

                        if job["status"] != current_status:
                            job["status"] = current_status
                            save_submitted_job(JOB_TRACKER_FILE, email, job)

                        if "CANCELLED" in current_status:
                            st.error("Job is cancelled. Contact support.")
                        elif current_status != "COMPLETED":
                            st.info(
                                "Job is still in progress. Check back later."
                            )
                        else:
                            st.success(
                                "Job is completed and ready for file retrieval."
                            )

        st.success("Job statuses checked successfully.")
    else:
        st.info("No jobs found for your email.")


def render_retrieve_output_tab(email):
    """
    Handles retrieving job output files.
    """
    st.header("Retrieve Job Output")
    job_id = st.text_input("Enter Job ID", key="job_id")

    jobs = load_submitted_jobs(JOB_TRACKER_FILE)
    job_name = None

    if job_id:
        for job in jobs.get(email, []):
            if job["job_id"] == job_id:
                job_name = job["job_name"]
                break

        if job_name:
            st.success(f"Job Name: {job_name}")
        else:
            st.warning("Job ID not found.")

    local_output_dir = "./output_files"

    if st.button("Retrieve Files", key="retrieve_files"):
        if not job_id or not job_name:
            st.warning("Please enter a valid Job ID.")
        else:
            server = RemoteServer(REMOTE_HOST, REMOTE_USER)
            job_status = JobStatus(server)

            status, error = job_status.check_status(job_id)
            if error:
                st.error(f"Error checking job status: {error}")
            elif status == "COMPLETED":
                job_folder_name = job_name.lower()
                job_folder_path = f"{local_output_dir}/{job_folder_name}"

                success, error = job_status.retrieve_files(
                    job_name, job_id, local_output_dir
                )
                if success:
                    st.success(
                        f"Files retrieved successfully and saved in {job_folder_path}"
                    )

                    zip_buffer = zip_job_folder(job_folder_path)
                    st.download_button(
                        label="Download Zipped Output Folder",
                        data=zip_buffer,
                        file_name=f"{job_folder_name}.zip",
                        mime="application/zip",
                    )
                else:
                    st.error(f"Error retrieving files: {error}")
            else:
                st.warning(
                    f"Job is not completed. Current status: {status}. Try later."
                )
