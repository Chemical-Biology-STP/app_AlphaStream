import os
import subprocess


class JobStatus:
    """Handles job status checking using SLURM's sacct command."""

    def __init__(self, server):
        self.server = server

    def check_status(self, job_id):
        """
        Check the real-time status of a job on the remote machine using SLURM.

        Args:
            job_id (str): The SLURM job ID.

        Returns:
            tuple: (status, error_message) where status is a string like "RUNNING" or "COMPLETED".
        """
        if not job_id or job_id == "Unknown":
            return None, "Invalid job ID."

        command = f"sacct -j {job_id} --format=JobID,State -n -P"
        output, error = self.server.run_command(command)

        if error:
            return None, error

        # Extract the job status from SLURM output
        lines = output.strip().split("\n")
        for line in lines:
            fields = line.split("|")
            if len(fields) > 1 and fields[0].strip() == job_id:
                return (
                    fields[1].strip(),
                    None,
                )  # Job status (e.g., PENDING, RUNNING, COMPLETED, FAILED)

        return None, "Job status not found."

    def retrieve_files(self, job_name, job_id, local_output_dir):
        """
        Retrieve output files from the remote machine to the local machine.

        Args:
            job_name (str): The job name.
            job_id (str): The SLURM job ID.
            local_output_dir (str): Local directory to store retrieved files.

        Returns:
            tuple: (success, error_message)
        """
        if not job_id or job_id == "Unknown":
            return False, "Invalid job ID."

        job_name_lower = job_name.lower()
        remote_output_path = f"{self.server.working_dir}/af_output/{job_name_lower}"  # ✅ Fix: Using self.server.working_dir
        local_output_path = os.path.join(local_output_dir, job_name_lower)

        # Ensure the local directory exists
        os.makedirs(local_output_path, exist_ok=True)

        # SCP command to copy the directory recursively
        scp_command = f"scp -r {self.server.user}@{self.server.host}:{remote_output_path} {local_output_path}"

        try:
            result = subprocess.run(
                scp_command, shell=True, capture_output=True, text=True
            )
            return (
                (True, None)
                if result.returncode == 0
                else (False, result.stderr.strip())
            )
        except Exception as e:
            return False, f"SCP error: {str(e)}"

    def get_all_active_jobs(self):
        """
        Retrieve all currently active jobs (PENDING/RUNNING) for the user.

        Returns:
            list: A list of dictionaries containing job details.
        """
        command = f"squeue -u {self.server.user} --format='%i|%T'"
        output, error = self.server.run_command(command)

        if error:
            return None, error

        active_jobs = []
        lines = output.strip().split("\n")[1:]  # Skip header
        for line in lines:
            job_id, state = line.split("|")
            active_jobs.append(
                {"job_id": job_id.strip(), "status": state.strip()}
            )

        return active_jobs, None
