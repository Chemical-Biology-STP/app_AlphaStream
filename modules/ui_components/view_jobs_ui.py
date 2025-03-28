import streamlit as st
from config import REMOTE_HOST, REMOTE_USER
from modules.job_status import JobStatus
from modules.remote_server import RemoteServer
from modules.utils import load_submitted_jobs, save_submitted_job


def render_view_jobs_tab(email):
    """
    Displays submitted jobs and provides options to refresh statuses.
    """
    st.header("View Submitted Jobs")

    jobs = load_submitted_jobs()
    user_jobs = jobs.get(email, [])

    if not user_jobs:
        st.info("No jobs found.")
        return

    server = RemoteServer(REMOTE_HOST, REMOTE_USER)
    job_status_checker = JobStatus(server)

    # Filter jobs that are still running or pending
    active_jobs = [
        job for job in user_jobs if job["status"] in ["Pending", "RUNNING"]
    ]

    # "Refresh All" button to check all running/pending jobs
    if active_jobs and st.button("🔄 Refresh All", key="refresh_all_jobs"):
        with st.spinner("Updating all job statuses... Please wait."):
            for job in active_jobs:
                job_id = job["job_id"]
                current_status, error = job_status_checker.check_status(job_id)

                if error:
                    st.error(
                        f"Error checking status for {job['job_name']}: {error}"
                    )
                elif current_status:
                    job["status"] = current_status
                    save_submitted_job(email, job)

        st.success(
            "All jobs updated successfully! Refresh the page to see changes."
        )

    # Display job details
    for job in user_jobs:
        job_id = job.get("job_id")
        job_name = job.get("job_name")
        job_status = job.get("status")

        with st.expander(f"**Job: {job_name} (Status: {job_status})**"):
            st.write(f"**Job ID:** {job_id}")
            st.write(f"**Current Status:** {job_status}")

            # Individual Refresh button for each job
            if job_status in ["Pending", "RUNNING"]:
                if st.button(f"🔄 Refresh", key=f"refresh_{job_id}"):
                    with st.spinner(f"Checking status of {job_name}..."):
                        current_status, error = (
                            job_status_checker.check_status(job_id)
                        )

                        if error:
                            st.error(f"Error checking status: {error}")
                        elif current_status:
                            job["status"] = current_status
                            save_submitted_job(email, job)
                            st.success(f"Updated status: {current_status}")

            # Job is completed
            if job_status == "COMPLETED":
                st.success("✅ Job completed. You can retrieve the output.")

            # Job is cancelled
            elif "CANCELLED" in job_status:
                st.error("❌ Job was cancelled.")

            # Job is still running
            elif job_status in ["Pending", "RUNNING"]:
                st.info("🔄 Job is still in progress. Check back later.")

    st.success("View updated statuses or refresh as needed.")
