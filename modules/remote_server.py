import subprocess

from config import REMOTE_WORKING_DIR


class RemoteServer:
    """Handles operations with the remote server using SSH and SCP."""

    def __init__(self, host, user):
        self.host = host
        self.user = user
        self.working_dir = (
            REMOTE_WORKING_DIR  # ✅ Stores remote working directory
        )

    def run_command(self, command):
        """
        Runs a command on the remote server via SSH.

        Args:
            command (str): The command to execute remotely.

        Returns:
            tuple: (output, error_message) where output is the command result.
        """
        try:
            ssh_command = f"ssh {self.user}@{self.host} '{command}'"
            result = subprocess.run(
                ssh_command,
                shell=True,
                capture_output=True,
                text=True,
                timeout=2400,
            )
            if result.returncode == 0:
                return result.stdout.strip(), None
            else:
                return None, result.stderr.strip()
        except subprocess.TimeoutExpired:
            return None, "SSH error: Command timed out"
        except Exception as e:
            return None, f"SSH error: {str(e)}"

    def transfer_file(self, local_path, remote_path):
        """
        Transfers a file to the remote server using SCP.

        Args:
            local_path (str): Path of the local file to send.
            remote_path (str): Destination path on the remote server.

        Returns:
            tuple: (success, error_message)
        """
        try:
            scp_command = (
                f"scp -r {local_path} {self.user}@{self.host}:{remote_path}"
            )
            result = subprocess.run(
                scp_command, shell=True, capture_output=True, text=True
            )
            return result.returncode == 0, (
                None if result.returncode == 0 else result.stderr.strip()
            )
        except Exception as e:
            return False, f"SCP error: {str(e)}"

    def check_file_exists(self, remote_path):
        """
        Checks if a file exists on the remote server.

        Args:
            remote_path (str): Path of the file on the remote server.

        Returns:
            bool: True if the file exists, False otherwise.
        """
        command = f"test -f {remote_path} && echo 'exists'"
        output, error = self.run_command(command)
        return output.strip() == "exists" if output else False

    def create_remote_directory(self, remote_path):
        """
        Creates a directory on the remote server if it does not already exist.

        Args:
            remote_path (str): Path of the directory to create.

        Returns:
            tuple: (success, error_message)
        """
        command = f"mkdir -p {remote_path}"
        output, error = self.run_command(command)
        return output is not None, error
