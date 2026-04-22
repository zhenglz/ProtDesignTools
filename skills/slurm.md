# SLURM Skill

## Overview
SLURM (Simple Linux Utility for Resource Management) is a job scheduler for Linux clusters. This document covers SLURM usage specific to this system, including available partitions, job submission, and script writing.

## System-Specific Partitions

### GPU Partitions
1. **4090**: NVIDIA RTX 4090 GPUs
   - High-performance consumer GPUs
   - 24GB VRAM per GPU
   - Best for most ML/AI workloads

2. **gpu_part**: General GPU partition
   - Mixed GPU types (may include A100, V100, etc.)
   - Default GPU partition

3. **p100**: NVIDIA P100 GPUs
   - Older data center GPUs
   - 16GB VRAM per GPU
   - Good for stable, well-tested workloads

### CPU Partitions
1. **AMD**: Modern AMD CPU nodes
   - Latest AMD EPYC processors
   - High core count, good for parallel CPU tasks

2. **old_CPU**: Older CPU nodes
   - Legacy Intel/AMD processors
   - Lower priority, longer wait times
   - Use for non-urgent batch jobs

## Basic SLURM Commands

### Job Submission
```bash
# Submit a job script
sbatch job_script.sh

# Submit with specific partition
sbatch -p 4090 job_script.sh

# Submit with job name
sbatch -J "protein_folding" job_script.sh
```

### Job Monitoring
```bash
# View all jobs
squeue

# View your jobs
squeue -u $USER

# View specific partition
squeue -p 4090

# View job details
scontrol show job <job_id>

# View job output while running
tail -f slurm-<job_id>.out
```

### Job Management
```bash
# Cancel a job
scancel <job_id>

# Cancel all your jobs
scancel -u $USER

# Hold a job
scontrol hold <job_id>

# Release a held job
scontrol release <job_id>
```

### Queue Information
```bash
# List all partitions
sinfo

# Detailed partition info
sinfo -o "%P %a %l %D %t %N %c %m %G"

# Check partition status
sinfo -p 4090,gpu_part,p100,AMD,old_CPU
```

## Writing SLURM Scripts

### Basic Template
```bash
#!/bin/bash
#SBATCH -J job_name              # Job name
#SBATCH -p partition_name        # Partition (4090, gpu_part, p100, AMD, old_CPU)
#SBATCH -N 1                     # Number of nodes
#SBATCH -n 1                     # Total number of tasks
#SBATCH --cpus-per-task=4        # CPUs per task
#SBATCH --mem-per-cpu=2GB        # Memory per CPU
#SBATCH --time=24:00:00          # Time limit (DD-HH:MM:SS)
#SBATCH -o slurm_output/%j.out   # Standard output
#SBATCH -e slurm_output/%j.err   # Standard error

# For GPU jobs
#SBATCH --gres=gpu:1             # Request 1 GPU
#SBATCH --gpus-per-task=1        # GPUs per task

echo "Starting job $SLURM_JOB_ID"
echo "Running on node: $SLURMD_NODENAME"
echo "Allocated GPUs: $CUDA_VISIBLE_DEVICES"

# Your commands here
module load compiler/gcc/11.2.1
module load compiler/cuda/12.2

python my_script.py

echo "Job completed"
```

### GPU-Specific Scripts

#### RTX 4090 Partition
```bash
#!/bin/bash
#SBATCH -J gpu_job_4090
#SBATCH -p 4090
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4GB
#SBATCH --gres=gpu:1
#SBATCH --time=48:00:00
#SBATCH -o %j.out
#SBATCH -e %j.err

# Load modules
module purge
module load compiler/gcc/11.2.1
module load compiler/cuda/12.4  # CUDA 12.4 for RTX 4090

# Set GPU environment
export CUDA_VISIBLE_DEVICES=0
export TF_FORCE_GPU_ALLOW_GROWTH=true

# Run GPU-intensive task
python deep_learning_model.py --batch_size 32 --epochs 100
```

#### General GPU Partition
```bash
#!/bin/bash
#SBATCH -J gpu_job_general
#SBATCH -p gpu_part
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00

module load compiler/gcc/11.2.1
module load compiler/cuda/12.2

# Detect GPU type
nvidia-smi --query-gpu=name --format=csv,noheader

# Run task
python protein_folding.py --model chai1
```

### CPU-Specific Scripts

#### AMD Partition (Modern CPUs)
```bash
#!/bin/bash
#SBATCH -J cpu_job_amd
#SBATCH -p AMD
#SBATCH -N 1
#SBATCH -n 32                     # Use many cores for parallel tasks
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=72:00:00

module load compiler/gcc/11.2.1

# Set OpenMP threads
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Run parallel CPU task
python molecular_dynamics.py --threads $SLURM_CPUS_PER_TASK
```

#### Old CPU Partition
```bash
#!/bin/bash
#SBATCH -J batch_job
#SBATCH -p old_CPU
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=16GB
#SBATCH --time=168:00:00          # 7 days for long-running batch jobs

# Non-urgent batch processing
python batch_processor.py --input data/ --output results/
```

## Using the Submission Script

The system provides a convenient submission script at `/data_test/home/lzzheng/bin/submit_slurm_gpu.sh`:

### Script Usage
```bash
# Basic usage
./submit_slurm_gpu.sh "<command>" <ncpus> <partition>

# Examples
./submit_slurm_gpu.sh "python train_model.py" 4 4090
./submit_slurm_gpu.sh "cd /path && ./run.sh" 8 gpu_part
./submit_slurm_gpu.sh "bash pipeline.sh" 2 p100
```

### Script Source
```bash
#!/bin/bash
# submit_slurm_gpu.sh

if [ $# -lt 2 ]; then
  echo "usage: submit_slurm.sh CMD NCPUS Queue"
  echo "example: submit_slurm_gpu.sh \"ls -rlt\" 4 4090"
  exit 0;
fi

cmd=${1}
ncpus=${2}

if [ $# -eq 3 ]; then
  queue="$3"
else
  queue="gpu_part"  # Default partition
fi

# Generate random script name
rndstr=$(head /dev/urandom | tr -dc A-Za-z0-9 | head -c10)

mkdir -pv slurm_output

cat <<EOT >> run_${rndstr}.sh
#!/bin/bash
#SBATCH -J general
#SBATCH -p ${queue}
#SBATCH -N 1
#SBATCH --time=96:00:00          # 4 days maximum
#SBATCH --mem-per-cpu=2GB
#SBATCH --ntasks-per-node=${ncpus}
#SBATCH -n ${ncpus}
#SBATCH --gres=gpu:1             # Always requests 1 GPU
#SBATCH --get-user-env
#SBATCH -e slurm_output/%j.err
#SBATCH -o slurm_output/%j.out

module load compiler/gcc/11.2.1
module list

date
time ${cmd}
date

echo "COMPLETE ..."
EOT

sbatch run_${rndstr}.sh
rm -rfv run_${rndstr}.sh
```

## Environment Modules

### Available Modules
```bash
# List available modules
module avail

# Load modules
module load compiler/gcc/11.2.1
module load compiler/cuda/12.2
module load compiler/cuda/12.4  # For RTX 4090

# List loaded modules
module list

# Unload modules
module unload compiler/cuda

# Purge all modules
module purge
```

### Common Module Combinations
```bash
# For GPU jobs (general)
module load compiler/gcc/11.2.1
module load compiler/cuda/12.2

# For RTX 4090 jobs
module load compiler/gcc/11.2.1
module load compiler/cuda/12.4

# For CPU-only jobs
module load compiler/gcc/11.2.1
```

## Resource Management

### Estimating Resources
```bash
# Check job efficiency after completion
seff <job_id>

# Example output:
# Job ID: 123456
# Cluster: cluster_name
# User/Group: user/group
# State: COMPLETED (exit code 0)
# Nodes: 1
# Cores per node: 8
# CPU Utilized: 2-04:12:34
# CPU Efficiency: 85.65% of 2-06:40:00 core-walltime
# Job Wall-clock time: 7-12:00:00
# Memory Utilized: 15.67 GB
# Memory Efficiency: 48.97% of 32.00 GB
```

### Resource Limits
```bash
# Check your limits
sacctmgr show user $USER

# Check partition limits
scontrol show partition <partition_name>
```

## Best Practices

### 1. Request Appropriate Resources
```bash
# Don't over-request
# BAD: Requesting 128GB for a 1GB job
# GOOD: Estimate based on previous runs

# Use seff to optimize future requests
seff <completed_job_id>
```

### 2. Use Output Directories
```bash
# Create output directory
mkdir -p slurm_output

# Redirect output
#SBATCH -o slurm_output/%j.out
#SBATCH -e slurm_output/%j.err
```

### 3. Include Module Loads
```bash
# Always load required modules
module load compiler/gcc/11.2.1
module load compiler/cuda/12.2
```

### 4. Set Time Limits Appropriately
```bash
# Short test jobs: 1-2 hours
#SBATCH --time=02:00:00

# Medium jobs: 1 day
#SBATCH --time=24:00:00

# Long jobs: 3-7 days
#SBATCH --time=168:00:00
```

### 5. Use Job Arrays for Similar Tasks
```bash
#!/bin/bash
#SBATCH -J array_job
#SBATCH -p AMD
#SBATCH -a 1-100                 # 100 array tasks
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2GB

# Each array task gets a different input
input_file="data/input_${SLURM_ARRAY_TASK_ID}.fasta"
output_file="results/output_${SLURM_ARRAY_TASK_ID}.txt"

python process.py --input "$input_file" --output "$output_file"
```

## Troubleshooting

### Common Issues
```bash
# Job pending forever
# Check: sinfo -p <partition> (partition may be down)
# Solution: Try different partition or contact admin

# Job failing immediately
# Check: tail -f slurm-<job_id>.err
# Common causes: module not loaded, path error

# Out of memory
# Check: seff <job_id> for memory usage
# Solution: Increase --mem or --mem-per-cpu

# GPU errors
# Check: nvidia-smi in job script
# Solution: Load correct CUDA version
```

### Debug Commands
```bash
# Test SLURM connectivity
sinfo

# Check your active jobs
squeue -u $USER -o "%.18i %.9P %.8j %.8u %.2t %.10M %.6D %.4C %.6m %R"

# Check job history
sacct -u $USER --format=JobID,JobName,Partition,AllocCPUS,State,ExitCode

# Check specific job
scontrol show job <job_id>
```

## Example Workflows

### Protein Folding with Chai-1
```bash
#!/bin/bash
#SBATCH -J chai1_folding
#SBATCH -p 4090
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4GB
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00

module load compiler/gcc/11.2.1
module load compiler/cuda/12.4

cd /data_test/share/pub_tools/chai-lab
./run.sh input.fasta ./output_chai1
```

### Molecular Dynamics with OpenMM
```bash
#!/bin/bash
#SBATCH -J openmm_simulation
#SBATCH -p gpu_part
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4GB
#SBATCH --gres=gpu:1
#SBATCH --time=72:00:00

module load compiler/gcc/11.2.1
module load compiler/cuda/12.2

python run_simulation.py --input protein.pdb --length 100 --platform CUDA
```

### Batch Sequence Processing
```bash
#!/bin/bash
#SBATCH -J batch_processing
#SBATCH -p AMD
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=48:00:00
#SBATCH -a 1-50

module load compiler/gcc/11.2.1

input_seq="sequences/seq_${SLURM_ARRAY_TASK_ID}.fasta"
output_dir="results/run_${SLURM_ARRAY_TASK_ID}"

python process_sequence.py --input "$input_seq" --output "$output_dir"
```

## References
- SLURM documentation: `man sbatch`, `man scontrol`
- System partitions: 4090, gpu_part, p100, AMD, old_CPU
- Submission script: `/data_test/home/lzzheng/bin/submit_slurm_gpu.sh`
- Module system: `module avail`, `module load`