import json
import os
import zipfile
from io import BytesIO

from config import JOB_TRACKER_FILE


def load_submitted_jobs():
    """
    Loads the submitted jobs from the job tracker file.
    Returns an empty dictionary if the file does not exist.
    """
    if os.path.exists(JOB_TRACKER_FILE):
        with open(JOB_TRACKER_FILE, "r") as f:
            return json.load(f)
    return {}


def save_submitted_job(email, job_data):
    """
    Saves job details to the tracker file. Updates existing jobs if needed.
    """
    jobs = load_submitted_jobs()

    if email not in jobs:
        jobs[email] = []

    # Update existing job if it already exists
    for job in jobs[email]:
        if job["job_name"] == job_data["job_name"]:
            job.update(job_data)
            break
    else:
        jobs[email].append(job_data)

    with open(JOB_TRACKER_FILE, "w") as f:
        json.dump(jobs, f, indent=4)


def remove_empty_fields(data):
    """
    Recursively removes keys with empty, null, or default values from a dictionary or list.
    """
    if isinstance(data, dict):
        return {
            k: remove_empty_fields(v)
            for k, v in data.items()
            if v not in [None, "", [], {}, 0]
        }
    elif isinstance(data, list):
        return [
            remove_empty_fields(v)
            for v in data
            if v not in [None, "", [], {}, 0]
        ]
    return data


def zip_job_folder(job_folder_path):
    """
    Compresses a specific job folder into a zip file and returns the zip file as a BytesIO object.
    """
    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zipf:
        for root, _, files in os.walk(job_folder_path):
            for file in files:
                file_path = os.path.join(root, file)
                arcname = os.path.relpath(file_path, job_folder_path)
                zipf.write(file_path, arcname)
    zip_buffer.seek(0)
    return zip_buffer


def validate_json_file(json_content):
    """
    Validates a JSON object to ensure proper formatting.

    Returns:
        (bool, str): True if valid, otherwise False with an error message.
    """
    try:
        json.loads(
            json.dumps(json_content)
        )  # Ensures JSON serialization is valid
        return True, ""
    except json.JSONDecodeError as e:
        return False, f"Invalid JSON format: {str(e)}"


def generate_unique_job_id():
    """
    Generates a unique job ID for tracking.
    """
    import uuid

    return str(uuid.uuid4())


def get_job_status(email, job_name):
    """
    Retrieves the status of a submitted job.
    """
    jobs = load_submitted_jobs()
    if email in jobs:
        for job in jobs[email]:
            if job["job_name"] == job_name:
                return job.get("status", "Unknown")
    return "Job not found"
