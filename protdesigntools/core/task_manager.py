
import os
import subprocess
import logging
import time
import json
from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Any, Union
from enum import Enum

logger = logging.getLogger(__name__)

class JobStatus(Enum):
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    COMPLETED = "COMPLETED"
    FAILED = "FAILED"
    CANCELLED = "CANCELLED"
    TIMEOUT = "TIMEOUT"
    UNKNOWN = "UNKNOWN"

class TaskManager(ABC):
    """Base class for task management (Local or Slurm)"""
    
    def __init__(self, max_jobs: int = 10, work_dir: str = "./work_dir"):
        self.max_jobs = max_jobs
        self.work_dir = work_dir
        os.makedirs(work_dir, exist_ok=True)
        self.jobs: Dict[str, Dict[str, Any]] = {}

    @abstractmethod
    def submit(self, command: str, job_name: str, **kwargs) -> str:
        """Submit a job and return job_id"""
        pass

    @abstractmethod
    def get_status(self, job_id: str) -> JobStatus:
        """Get status of a specific job"""
        pass

    @abstractmethod
    def cancel(self, job_id: str) -> bool:
        """Cancel a specific job"""
        pass

    def wait_for_jobs(self, job_ids: List[str], poll_interval: int = 10) -> None:
        """Wait for a list of jobs to complete"""
        while True:
            all_done = True
            for job_id in job_ids:
                status = self.get_status(job_id)
                if status in (JobStatus.PENDING, JobStatus.RUNNING):
                    all_done = False
                    break
            if all_done:
                break
            time.sleep(poll_interval)

class LocalTaskManager(TaskManager):
    """Local task manager using subprocess"""
    
    def __init__(self, max_jobs: int = 4, work_dir: str = "./work_dir", **kwargs):
        super().__init__(max_jobs, work_dir)
        self.processes: Dict[str, subprocess.Popen] = {}

    def submit(self, command: str, job_name: str, **kwargs) -> str:
        # Simple local execution. In a real scenario, we might want a queue system.
        # For now, we just run it and store the process object.
        log_file = os.path.join(self.work_dir, f"{job_name}.log")
        
        # We need to keep the file object open while the process is running,
        # otherwise stdout/stderr writes will fail.
        f = open(log_file, "w")
        proc = subprocess.Popen(
            command, 
            shell=True, 
            stdout=f, 
            stderr=subprocess.STDOUT,
            cwd=self.work_dir
        )
        
        # Store the file object so we can close it later when process finishes
        proc._log_file_obj = f
        
        job_id = str(proc.pid)
        self.processes[job_id] = proc
        self.jobs[job_id] = {
            "name": job_name,
            "command": command,
            "status": JobStatus.RUNNING,
            "log": log_file
        }
        return job_id

    def get_status(self, job_id: str) -> JobStatus:
        if job_id not in self.processes:
            return JobStatus.UNKNOWN
        
        proc = self.processes[job_id]
        ret = proc.poll()
        if ret is None:
            return JobStatus.RUNNING
        
        if hasattr(proc, "_log_file_obj") and not proc._log_file_obj.closed:
            proc._log_file_obj.close()
            
        if ret == 0:
            return JobStatus.COMPLETED
        else:
            return JobStatus.FAILED

    def cancel(self, job_id: str) -> bool:
        if job_id in self.processes:
            proc = self.processes[job_id]
            proc.terminate()
            if hasattr(proc, "_log_file_obj"):
                proc._log_file_obj.close()
            return True
        return False
        
    def wait_for_jobs(self, job_ids: List[str], poll_interval: int = 10) -> None:
        """Wait for a list of jobs to complete"""
        for job_id in job_ids:
            if job_id in self.processes:
                proc = self.processes[job_id]
                # wait() can return returncode, we wait for the process to finish
                proc.wait()
                
                # close the file handler
                if hasattr(proc, "_log_file_obj"):
                    proc._log_file_obj.close()
                    
                # Check if it succeeded
                if proc.returncode != 0:
                    logger.warning(f"Local job {job_id} exited with code {proc.returncode}")

class SlurmTaskManager(TaskManager):
    """Slurm task manager"""
    
    def __init__(self, max_jobs: int = 100, work_dir: str = "./work_dir", partition: str = "AMD", **kwargs):
        super().__init__(max_jobs, work_dir)
        self.partition = partition

    def submit(self, command: str, job_name: str, **kwargs) -> str:
        partition = kwargs.get("partition", self.partition)
        nodes = kwargs.get("nodes", 1)
        ntasks = kwargs.get("ntasks", 1)
        cpus_per_task = kwargs.get("cpus_per_task", 4)
        gres = kwargs.get("gres", "") # e.g. gpu:1
        
        slurm_script = os.path.join(self.work_dir, f"{job_name}.sh")
        with open(slurm_script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(f"#SBATCH --job-name={job_name}\n")
            f.write(f"#SBATCH --partition={partition}\n")
            f.write(f"#SBATCH --nodes={nodes}\n")
            f.write(f"#SBATCH --ntasks={ntasks}\n")
            f.write(f"#SBATCH --cpus-per-task={cpus_per_task}\n")
            if gres:
                f.write(f"#SBATCH --gres={gres}\n")
            f.write(f"#SBATCH --output={job_name}.out\n")
            f.write(f"#SBATCH --error={job_name}.err\n")
            f.write(f"\n{command}\n")
        
        result = subprocess.run(f"sbatch {slurm_script}", shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            # Example output: "Submitted batch job 12345"
            job_id = result.stdout.strip().split()[-1]
            self.jobs[job_id] = {
                "name": job_name,
                "command": command,
                "status": JobStatus.PENDING
            }
            return job_id
        else:
            logger.error(f"Slurm submission failed: {result.stderr}")
            raise RuntimeError(f"Failed to submit slurm job: {result.stderr}")

    def get_status(self, job_id: str) -> JobStatus:
        result = subprocess.run(f"squeue -j {job_id} -h -o %T", shell=True, capture_output=True, text=True)
        status_str = result.stdout.strip()
        
        if not status_str:
            # Check sacct for finished jobs
            result = subprocess.run(f"sacct -j {job_id} -h -o State", shell=True, capture_output=True, text=True)
            status_str = result.stdout.strip().split('\n')[0].split()[0] if result.stdout.strip() else ""
            
        if "PENDING" in status_str: return JobStatus.PENDING
        if "RUNNING" in status_str: return JobStatus.RUNNING
        if "COMPLETED" in status_str: return JobStatus.COMPLETED
        if "FAILED" in status_str: return JobStatus.FAILED
        if "CANCELLED" in status_str: return JobStatus.CANCELLED
        if "TIMEOUT" in status_str: return JobStatus.TIMEOUT
        return JobStatus.UNKNOWN

    def cancel(self, job_id: str) -> bool:
        result = subprocess.run(f"scancel {job_id}", shell=True)
        return result.returncode == 0

def get_manager(mode: str = "local", **kwargs) -> TaskManager:
    if mode.lower() == "slurm":
        return SlurmTaskManager(**kwargs)
    else:
        return LocalTaskManager(**kwargs)
