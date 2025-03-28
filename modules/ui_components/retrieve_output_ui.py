import os

import streamlit as st
from config import REMOTE_HOST, REMOTE_USER
from modules.job_status import JobStatus
from modules.remote_server import RemoteServer
from modules.utils import load_submitted_jobs, zip_job_folder


def render_retrieve_output_tab(email):
    """
    Handles retrieving job output files and provides a download link.
    """
    st.header("Retrieve Job Output")

    # User inputs job ID
    job_id = st.text_input("Enter Job ID", key="job_id")

    # Load submitted jobs
    jobs = load_submitted_jobs()
    job_name = None

    if job_id:
        # Find job associated with user and entered Job ID
        for job in jobs.get(email, []):
            if job["job_id"] == job_id:
                job_name = job["job_name"]
                break

        if not job_name:
            st.warning("⚠️ Job ID not found in your submitted jobs.")
            return  # ✅ Fix: Stop execution if job_name is None

    # ✅ Fix: Ensure job_name is defined before using `.lower()`
    job_folder_name = job_name.lower() if job_name else None
    local_output_dir = "./output_files"
    job_folder_path = (
        os.path.join(local_output_dir, job_folder_name)
        if job_folder_name
        else None
    )

    # If job_name is None, exit the function
    if job_folder_name is None or job_folder_path is None:
        return

    # Check if files are already downloaded
    if os.path.exists(job_folder_path) and os.listdir(job_folder_path):
        # Zip the existing job folder for downloading
        zip_buffer = zip_job_folder(job_folder_path)

        # Provide download button immediately
        st.download_button(
            label="⬇️ Download Job Output (ZIP)",
            data=zip_buffer,
            file_name=f"{job_folder_name}.zip",
            mime="application/zip",
        )
    else:
        # If files are not found locally, provide the retrieve button
        if st.button("📂 Retrieve Files", key="retrieve_files"):
            if not job_id:
                st.warning("⚠️ Please enter a Job ID.")
                return

            server = RemoteServer(REMOTE_HOST, REMOTE_USER)
            job_status = JobStatus(server)

            # Check job status before retrieving files
            status, error = job_status.check_status(job_id)
            if error:
                st.error(f"❌ Error checking job status: {error}")
                return
            elif status != "COMPLETED":
                st.warning(
                    f"⏳ Job is not completed yet. Current status: {status}. Try again later."
                )
                return

            # Retrieve job output files
            success, error = job_status.retrieve_files(
                job_name, job_id, local_output_dir
            )
            if success:
                # Zip the newly retrieved job folder for downloading
                zip_buffer = zip_job_folder(job_folder_path)

                # Provide download button
                st.download_button(
                    label="⬇️ Download Job Output (ZIP)",
                    data=zip_buffer,
                    file_name=f"{job_folder_name}.zip",
                    mime="application/zip",
                )
            else:
                st.error(f"❌ Error retrieving files: {error}")
