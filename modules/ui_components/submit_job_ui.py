import json

import streamlit as st
from config import REMOTE_HOST, REMOTE_USER
from modules.input_validation import is_valid_string
from modules.job_submission import JobSubmission
from modules.ligand_handler import collect_ligands
from modules.nucleic_acid_handler import collect_nucleic_acid_sequences
from modules.protein_handler import collect_protein_sequences
from modules.remote_server import RemoteServer
from modules.utils import remove_empty_fields, save_submitted_job

# ✅ Load CCD.json (Modify path if necessary)
#with open("CCD.json", "r") as f:
with open("app_AlphaStream/CCD.json", "r") as f:
    CCD_DATA = json.load(f)

# ✅ Allowed characters for DNA & RNA validation
VALID_DNA_CHARS = set("ACGT")
VALID_RNA_CHARS = set("ACGU")


def render_submit_job_tab(email):
    """
    Renders the UI for job submission, including Protein, Ligand, DNA & RNA inputs.
    """
    st.header("Submit Job")
    job_name = st.text_input("Job Name")

    # Validate Job Name
    valid_job_name, error_message = is_valid_string(job_name)
    if not valid_job_name:
        st.error(
            "❌ Invalid Job Name. Only alphabets, numbers, dashes, and underscores are allowed."
        )
    if not job_name:
        st.error("❌ Job Name cannot be empty.")

    model_seeds = st.text_area("Model Seeds (comma-separated)", value="1")
    unique_ids = set()

    # Collect Protein Sequences
    protein_sequences, proteins_valid = collect_protein_sequences(unique_ids)

    # Collect Ligands
    ligands, ligands_valid = collect_ligands(unique_ids)

    # Collect DNA & RNA Sequences
    dna_sequences, dna_valid = collect_nucleic_acid_sequences(
        "DNA", VALID_DNA_CHARS, unique_ids
    )
    rna_sequences, rna_valid = collect_nucleic_acid_sequences(
        "RNA", VALID_RNA_CHARS, unique_ids
    )

    # Final Input Validation
    all_inputs_valid = (
        valid_job_name
        and job_name
        and proteins_valid
        and ligands_valid
        and dna_valid
        and rna_valid
    )

    # Construct JSON payload
    json_content = {
        "name": job_name,
        "modelSeeds": [
            int(seed.strip()) for seed in model_seeds.split(",") if seed.strip()
        ],
        "sequences": protein_sequences + ligands + dna_sequences + rna_sequences,
        "dialect": "alphafold3",
        "version": 1,
    }
    json_content = remove_empty_fields(json_content)

    st.subheader("📜 Generated JSON")
    st.json(json_content)

    st.download_button(
        "⬇️ Download JSON File",
        json.dumps(json_content, indent=4),
        f"{job_name}.json",
        "application/json",
        disabled=not all_inputs_valid,
    )
    if st.button("🚀 Submit Job", key="submit_job", disabled=not all_inputs_valid):
        handle_job_submission(email, job_name, json_content)


def handle_job_submission(email, job_name, json_content):
    """Handles background job submission."""
    save_submitted_job(
        email, {"job_id": "Pending", "job_name": job_name, "status": "Pending"}
    )
    JobSubmission(RemoteServer(REMOTE_HOST, REMOTE_USER)).submit(
        json_content, email, job_name
    )
    st.success("✅ Your job is being submitted in the background.")
