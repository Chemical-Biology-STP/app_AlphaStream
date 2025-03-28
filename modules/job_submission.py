import json
import os
import re
from uuid import uuid4

from config import REMOTE_WORKING_DIR
from modules.remote_server import RemoteServer
from modules.utils import generate_unique_job_id, save_submitted_job


class JobSubmission:
    """Handles job submission logic, including remote file transfer and SLURM job submission."""

    def __init__(self, server: RemoteServer):
        self.server = server

    def submit(self, json_content, email, job_name):
        """
        Submits a job by transferring input files to the remote server and executing a SLURM job.

        Args:
            json_content (dict): The job's input JSON.
            email (str): User's email address.
            job_name (str): The name of the job.

        Returns:
            tuple: (job_id, error_message) where job_id is the SLURM job ID if successful, otherwise None.
        """
        # Generate a local JSON file
        local_json_file = f"{job_name}.json"
        with open(local_json_file, "w") as f:
            json.dump(json_content, f, indent=4)

        # Define remote paths
        remote_json_path = os.path.join(
            REMOTE_WORKING_DIR, "af_input", local_json_file
        )

        # Transfer JSON file to remote server
        transfer_status, transfer_error = self.server.transfer_file(
            local_json_file, remote_json_path
        )
        if not transfer_status:
            return None, f"File transfer failed: {transfer_error}"

        # Generate SLURM script
        job_id = generate_unique_job_id()
        bash_script_filename = f"af3_{job_id}.sh"
        bash_script_content = self._generate_bash_script(
            email, job_name, remote_json_path
        )

        # Save the SLURM script locally
        with open(bash_script_filename, "w") as script_file:
            script_file.write(bash_script_content)

        # Transfer SLURM script to remote server
        remote_script_path = os.path.join(
            REMOTE_WORKING_DIR, "af_input", bash_script_filename
        )
        transfer_status, transfer_error = self.server.transfer_file(
            bash_script_filename, remote_script_path
        )
        if not transfer_status:
            return None, f"SLURM script transfer failed: {transfer_error}"

        # Submit SLURM job
        result_output, result_error = self.server.run_command(
            f"sbatch {remote_script_path}"
        )

        # Clean up local temporary files
        os.remove(local_json_file)
        os.remove(bash_script_filename)

        if result_output:
            # Extract SLURM job ID from submission response
            match = re.search(r"Submitted batch job (\d+)", result_output)
            slurm_job_id = match.group(1) if match else "Unknown"

            # Store the submitted job in the tracker
            job_data = {
                "job_id": slurm_job_id,
                "job_name": job_name,
                "status": "Pending",
            }
            save_submitted_job(email, job_data)

            return slurm_job_id, None  # Successfully submitted, return job ID
        else:
            job_data = {
                "job_id": "Failed",
                "job_name": job_name,
                "status": "Failed",
            }
            save_submitted_job(email, job_data)
            return None, f"Job submission failed: {result_error}"

    @staticmethod
    def _generate_bash_script(email, job_name, json_path):
        """
        Generates a SLURM batch script for submitting the job.

        Args:
            email (str): User's email.
            job_name (str): The job name.
            json_path (str): Path to the JSON input file.

        Returns:
            str: SLURM batch script content.
        """
        return f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --time=1-00:00:00
#SBATCH --partition=gh100
#SBATCH --gres=gpu:1
#SBATCH --mem=1000G
#SBATCH --reservation=h100
#SBATCH --output={REMOTE_WORKING_DIR}/af_output/{job_name}_%j.out
#SBATCH --error={REMOTE_WORKING_DIR}/af_output/{job_name}_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user={email}

ml Singularity
AF="{REMOTE_WORKING_DIR}"

singularity exec \\
     --nv \\
     --bind $AF/af_input:/root/af_input \\
     --bind {REMOTE_WORKING_DIR}/af_output:/root/af_output \\
     --bind $AF/models:/root/models \\
     --bind $AF/public_databases:/root/public_databases \\
     /nemo/stp/chemicalbiology/home/users/yipy/Documents/alphafold3.sif \\
     python alphafold3/run_alphafold.py \\
     --json_path=/root/af_input/{os.path.basename(json_path)} \\
     --model_dir=/root/models \\
     --db_dir=/root/public_databases \\
     --output_dir=/root/af_output \\
     --flash_attention_implementation=xla

# Convert SLURM_JOB_NAME to lowercase
lowercase_job_name=$(echo "$SLURM_JOB_NAME" | tr '[:upper:]' '[:lower:]')

# Define the link
link="/nemo/stp/chemicalbiology/home/users/yipy/Documents/af3/af_output/$lowercase_job_name"

# Compose the email
TO="{email}"
SUBJECT="Job Completed: $lowercase_job_name"

(
  echo "To: $TO"
  echo "Subject: $SUBJECT"
  echo "Content-Type: text/html"
  echo
  echo "<html>"
  echo "<body>"
  echo "<p>Job completed. Please download the files <a href='file://$link'>here</a>.</p>"
  echo "</body>"
  echo "</html>"
) | sendmail -t
"""
