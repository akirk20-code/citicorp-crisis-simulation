# GMU Hopper HPC Setup Guide for OpenFOAM CFD

**Cluster:** GMU Hopper (https://wiki.orc.gmu.edu/mkdocs/About_Hopper/)
**Purpose:** Large-scale CFD runs (8-30M cells), parametric studies, LES simulations
**Target Hardware:** AMD EPYC 128-core nodes, NVIDIA A100 GPUs (optional)

---

## Table of Contents
1. [Account Setup & Access](#1-account-setup--access)
2. [Environment Configuration](#2-environment-configuration)
3. [Case Transfer & Preparation](#3-case-transfer--preparation)
4. [Slurm Job Scripts](#4-slurm-job-scripts)
5. [Monitoring & Troubleshooting](#5-monitoring--troubleshooting)
6. [Data Transfer Back to Windows](#6-data-transfer-back-to-windows)
7. [Cost Estimation & Allocation](#7-cost-estimation--allocation)

---

## 1. Account Setup & Access

### 1.1 Hopper Account Request
1. Visit: https://wiki.orc.gmu.edu/mkdocs/Applying_for_an_Account/
2. Faculty sponsor: [Your advisor's name]
3. Request compute allocation: **10,000 SUs (Service Units)** for initial testing
   - 1 SU = 1 core-hour
   - Example: 128-core job × 10 hours = 1,280 SUs

### 1.2 SSH Access from Windows

**Option A: WSL2 Ubuntu (Recommended)**
```bash
# From WSL2 terminal
ssh YOUR_GMU_ID@hopper.orc.gmu.edu

# Setup SSH key for password-less login
ssh-keygen -t ed25519 -C "your.email@gmu.edu"
ssh-copy-id YOUR_GMU_ID@hopper.orc.gmu.edu
```

**Option B: Windows PowerShell**
```powershell
ssh YOUR_GMU_ID@hopper.orc.gmu.edu
```

**Option C: VS Code Remote SSH (Best for editing)**
1. Install "Remote - SSH" extension
2. Add to `~/.ssh/config`:
```
Host hopper
    HostName hopper.orc.gmu.edu
    User YOUR_GMU_ID
    ForwardX11 yes
    ServerAliveInterval 60
```
3. Connect: Ctrl+Shift+P → "Remote-SSH: Connect to Host" → hopper

### 1.3 Initial Login & Storage Check
```bash
# Check home directory quota (50 GB limit)
du -sh ~
quota -s

# Scratch directory (no quota, but purged after 30 days)
cd /scratch/YOUR_GMU_ID
mkdir -p openfoam_cases
cd openfoam_cases
```

**Storage Strategy:**
- **Home (~):** Store job scripts, small configs, final results (<5 GB)
- **Scratch (/scratch/YOUR_GMU_ID):** Run simulations, large time-series data
- **Transfer off-cluster:** Results > 30 days old (download to your laptop/OneDrive)

---

## 2. Environment Configuration

### 2.1 Check Available OpenFOAM Modules
```bash
module avail openfoam

# Expected output:
# openfoam/v2306
# openfoam/v2312
# openfoam/10
# openfoam/11
```

**Recommendation:** Use `openfoam/v2312` (ESI release, closest to your WSL ESI v2512)

### 2.2 Load OpenFOAM Environment
```bash
# Add to ~/.bashrc for automatic loading
echo "module load openfoam/v2312" >> ~/.bashrc
source ~/.bashrc

# Verify
which simpleFoam
# Should output: /apps/OpenFOAM/OpenFOAM-v2312/platforms/.../simpleFoam
```

### 2.3 MPI Configuration
```bash
# Hopper uses Intel MPI (preferred) or OpenMPI
module load impi/2021.5  # Intel MPI

# Verify MPI
mpirun --version
```

---

## 3. Case Transfer & Preparation

### 3.1 Transfer Case from WSL to Hopper

**Method 1: rsync (Recommended — preserves permissions, resumes on failure)**
```bash
# From WSL2 terminal
cd ~/citicorp_cfd
rsync -avz --progress \
    --exclude='*.foam' \
    --exclude='processor*' \
    --exclude='postProcessing' \
    --exclude='log.*' \
    --exclude='*.log' \
    ./ YOUR_GMU_ID@hopper.orc.gmu.edu:/scratch/YOUR_GMU_ID/openfoam_cases/citicorp_cfd/

# Explanation:
# -a: archive mode (preserves timestamps, permissions)
# -v: verbose
# -z: compress during transfer
# --exclude: skip large unnecessary files
```

**Method 2: scp (Simple, no resume)**
```bash
cd ~
tar -czf citicorp_cfd.tar.gz citicorp_cfd/
scp citicorp_cfd.tar.gz YOUR_GMU_ID@hopper.orc.gmu.edu:/scratch/YOUR_GMU_ID/

# On Hopper:
ssh hopper
cd /scratch/YOUR_GMU_ID
tar -xzf citicorp_cfd.tar.gz
```

**Method 3: Git (Best for version control)**
```bash
# On WSL, push to GitHub
cd ~/citicorp_cfd
git add -A
git commit -m "Prepare for Hopper HPC run"
git push origin main

# On Hopper, clone
ssh hopper
cd /scratch/YOUR_GMU_ID/openfoam_cases
git clone https://github.com/YOUR_USERNAME/citicorp-crisis-simulation.git citicorp_cfd
cd citicorp_cfd/citicorp_cfd
```

### 3.2 Clean Case for Fresh Run
```bash
# On Hopper
cd /scratch/YOUR_GMU_ID/openfoam_cases/citicorp_cfd

# Remove old results
rm -rf processor* postProcessing *.foam log.* dynamicCode

# Verify mesh exists (if not, will need to run meshing first)
ls -lh constant/polyMesh/
# Should see: points, faces, owner, neighbour, boundary
```

---

## 4. Slurm Job Scripts

### 4.1 Basic Slurm Job Script (RANS, 3M cells, 64 cores)

Create `job_rans_3M.slurm`:
```bash
#!/bin/bash
#SBATCH --job-name=citicorp_rans_3M
#SBATCH --partition=normal         # normal, large-mem, or gpuq
#SBATCH --nodes=1                  # Single node (up to 128 cores)
#SBATCH --ntasks=64                # 64 MPI processes
#SBATCH --cpus-per-task=1          # 1 thread per process
#SBATCH --mem=128G                 # 128 GB RAM (2 GB per core)
#SBATCH --time=12:00:00            # 12 hour walltime
#SBATCH --output=slurm_%j.out      # stdout log (%j = job ID)
#SBATCH --error=slurm_%j.err       # stderr log
#SBATCH --mail-type=END,FAIL       # Email on completion or failure
#SBATCH --mail-user=YOUR_EMAIL@gmu.edu

# Load modules
module purge
module load openfoam/v2312
module load impi/2021.5

# Navigate to case directory
cd /scratch/$USER/openfoam_cases/citicorp_cfd

# Source OpenFOAM environment
source $FOAM_BASH

# Decompose mesh (if not already done)
if [ ! -d "processor0" ]; then
    echo "Decomposing mesh for 64 cores..."
    sed -i "s/numberOfSubdomains.*/numberOfSubdomains 64;/" system/decomposeParDict
    decomposePar > log.decomposePar 2>&1
fi

# Run solver
echo "Starting simpleFoam at $(date)"
mpirun -np $SLURM_NTASKS simpleFoam -parallel > log.simpleFoam 2>&1
echo "Finished simpleFoam at $(date)"

# Reconstruct final time only (save space)
echo "Reconstructing final timestep..."
reconstructPar -latestTime > log.reconstructPar 2>&1

# Calculate forces
echo "Calculating forces..."
simpleFoam -postProcess -func forces -latestTime >> log.forces 2>&1

# Cleanup: remove processor directories (save 10-20 GB)
# CAREFUL: Only if you don't need intermediate results
# rm -rf processor*

echo "Job completed successfully at $(date)"
```

**Submit:**
```bash
sbatch job_rans_3M.slurm
```

### 4.2 Mesh Generation Job (snappyHexMesh, 8M target)

Create `job_mesh_8M.slurm`:
```bash
#!/bin/bash
#SBATCH --job-name=citicorp_mesh_8M
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=16                # Fewer cores for meshing (less memory per core)
#SBATCH --mem=200G                 # High memory for snappyHexMesh
#SBATCH --time=6:00:00
#SBATCH --output=slurm_mesh_%j.out

module load openfoam/v2312 impi/2021.5
cd /scratch/$USER/openfoam_cases/citicorp_cfd_8M

source $FOAM_BASH

# Clean old mesh
rm -rf constant/polyMesh processor* 0/polyMesh

# Update blockMeshDict for finer base mesh (12m cells → 60×60×50 = 180k base)
# (Edit system/blockMeshDict manually before submitting)

# Update snappyHexMeshDict for higher refinement
# maxGlobalCells 10000000 (from 4M)
# maxLocalCells 5000000
# (Edit system/snappyHexMeshDict manually)

echo "Running blockMesh..."
blockMesh > log.blockMesh 2>&1

echo "Decomposing background mesh..."
decomposePar > log.decomposePar 2>&1

echo "Running snappyHexMesh (parallel, 16 cores)..."
mpirun -np 16 snappyHexMesh -parallel -overwrite > log.snappyHexMesh 2>&1

echo "Reconstructing mesh..."
reconstructParMesh -constant > log.reconstructParMesh 2>&1

echo "Checking mesh quality..."
checkMesh > log.checkMesh 2>&1

echo "Meshing complete. Check log.checkMesh for quality."
grep "Mesh OK" log.checkMesh && echo "SUCCESS" || echo "WARNING: Mesh issues detected"
```

**Submit:**
```bash
sbatch job_mesh_8M.slurm
```

### 4.3 LES Job (20M cells, 128 cores, GPU optional)

Create `job_les_20M.slurm`:
```bash
#!/bin/bash
#SBATCH --job-name=citicorp_les_20M
#SBATCH --partition=normal         # Or 'gpuq' for GPU nodes
#SBATCH --nodes=2                  # 2 nodes × 64 cores = 128 total
#SBATCH --ntasks=128
#SBATCH --mem-per-cpu=4G           # 4 GB × 128 = 512 GB total
#SBATCH --time=7-00:00:00          # 7 days walltime
#SBATCH --output=slurm_les_%j.out

module load openfoam/v2312 impi/2021.5
cd /scratch/$USER/openfoam_cases/citicorp_cfd_les

source $FOAM_BASH

# Decompose for 128 cores
if [ ! -d "processor0" ]; then
    sed -i "s/numberOfSubdomains.*/numberOfSubdomains 128;/" system/decomposeParDict
    decomposePar > log.decomposePar 2>&1
fi

# Run transient LES solver
echo "Starting pimpleFoam (LES) at $(date)"
mpirun -np 128 pimpleFoam -parallel > log.pimpleFoam 2>&1
echo "Finished at $(date)"

# Time-average last 50 seconds
echo "Time-averaging results..."
postProcess -func 'timeAverage(U,p)' -time '50:100' > log.timeAverage 2>&1

reconstructPar -latestTime > log.reconstructPar 2>&1
```

**Submit:**
```bash
sbatch job_les_20M.slurm
```

### 4.4 GPU-Accelerated Job (A100, 20M cells)

**Note:** Only use if Hopper has OpenFOAM with PETSc/AmgX GPU support (check with admins first)

Create `job_gpu_20M.slurm`:
```bash
#!/bin/bash
#SBATCH --job-name=citicorp_gpu
#SBATCH --partition=gpuq           # GPU partition
#SBATCH --nodes=1
#SBATCH --ntasks=8                 # Fewer CPU tasks, GPU does pressure solve
#SBATCH --gres=gpu:a100:1          # Request 1 A100 GPU
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm_gpu_%j.out

module load openfoam/v2312 cuda/12.1 impi/2021.5
cd /scratch/$USER/openfoam_cases/citicorp_cfd_gpu

source $FOAM_BASH

# Verify GPU visible
nvidia-smi

# Run with GPU-accelerated pressure solver
decomposePar > log.decomposePar 2>&1
mpirun -np 8 simpleFoam -parallel > log.simpleFoam 2>&1
```

---

## 5. Monitoring & Troubleshooting

### 5.1 Job Status
```bash
# Check your jobs
squeue -u $USER

# Output format:
# JOBID  PARTITION  NAME          USER   ST  TIME  NODES  NODELIST
# 12345  normal     citicorp_rans kirka  R   2:30  1      node042

# Job states:
# PD = Pending (waiting for resources)
# R  = Running
# CG = Completing
# CD = Completed
# F  = Failed
```

### 5.2 Cancel Job
```bash
scancel 12345    # Cancel job 12345
scancel -u $USER # Cancel all your jobs
```

### 5.3 Monitor Running Job

**Tail the log file:**
```bash
tail -f /scratch/$USER/openfoam_cases/citicorp_cfd/log.simpleFoam
```

**Check convergence (while job is running):**
```bash
grep "Time =" log.simpleFoam | tail -20
```

**Monitor memory usage:**
```bash
sstat -j 12345 --format=JobID,MaxRSS,AvgRSS
```

**Check CPU efficiency:**
```bash
seff 12345  # Run after job completes

# Example output:
# Job ID: 12345
# Cluster: hopper
# User/Group: kirka/users
# Cores: 64
# CPU Utilized: 400:00:00
# CPU Efficiency: 78.12% (should be >70%)
# Memory Utilized: 120.5 GB
# Memory Efficiency: 94.14%
```

### 5.4 Common Errors & Fixes

**Error: `mpirun: command not found`**
```bash
# Fix: Load MPI module
module load impi/2021.5
```

**Error: `simpleFoam: command not found`**
```bash
# Fix: Source OpenFOAM environment
source $FOAM_BASH
# Or check module loaded:
module list
```

**Error: `Cannot find case directory`**
```bash
# Fix: Use absolute path in Slurm script
cd /scratch/$USER/openfoam_cases/citicorp_cfd
# Not: cd ~/openfoam_cases/citicorp_cfd (~ may not expand in batch mode)
```

**Error: Job killed with "OUT_OF_MEMORY"**
```bash
# Fix: Request more memory
#SBATCH --mem=256G   # Increase from 128G

# Or reduce core count (more RAM per core)
#SBATCH --ntasks=32  # Instead of 64
```

**Error: `decomposeParDict: numberOfSubdomains 20 != ntasks 64`**
```bash
# Fix: Update decomposeParDict before submitting
sed -i "s/numberOfSubdomains.*/numberOfSubdomains 64;/" system/decomposeParDict
```

---

## 6. Data Transfer Back to Windows

### 6.1 Transfer Results to WSL

**Method 1: rsync (resume-able)**
```bash
# From WSL2 terminal
rsync -avz --progress \
    YOUR_GMU_ID@hopper.orc.gmu.edu:/scratch/YOUR_GMU_ID/openfoam_cases/citicorp_cfd/3000/ \
    ~/citicorp_cfd_results/3000/
```

**Method 2: scp (single file)**
```bash
scp YOUR_GMU_ID@hopper.orc.gmu.edu:/scratch/YOUR_GMU_ID/.../log.simpleFoam ~/
```

**Method 3: Compress first (for large datasets)**
```bash
# On Hopper:
cd /scratch/$USER/openfoam_cases/citicorp_cfd
tar -czf results_3000.tar.gz 3000/ postProcessing/ log.*

# From WSL:
scp hopper:/scratch/$USER/openfoam_cases/citicorp_cfd/results_3000.tar.gz ~/
cd ~
tar -xzf results_3000.tar.gz
```

### 6.2 Transfer to OneDrive (Windows)

**Option A: Copy from WSL to Windows**
```bash
# From WSL
cp -r ~/citicorp_cfd_results /mnt/c/Users/kirka/OneDrive/Documents/GMU\ PhD/OR\ 750\ Reliability\,\ Safety\,\ and\ Risk/Citi\ Sim/
```

**Option B: Use OneDrive sync**
1. Move results to OneDrive folder in Windows Explorer
2. Wait for OneDrive sync to cloud
3. Access from any device

---

## 7. Cost Estimation & Allocation

### 7.1 SU (Service Unit) Calculation
```
SUs = cores × hours

Examples:
- 64 cores × 10 hours = 640 SUs
- 128 cores × 168 hours (7 days) = 21,504 SUs
```

### 7.2 Allocation Requests

**Initial allocation (free tier):** 10,000 SUs
- ~10 runs of 64-core jobs @ 10 hours each

**Medium allocation:** 50,000-100,000 SUs
- Parametric study: 20 wind directions × 3 wind speeds × 64 cores × 10 hrs = 38,400 SUs

**Large allocation (proposal required):** 500,000+ SUs
- LES campaign: 10 cases × 128 cores × 168 hrs = 215,040 SUs

**Request more:** https://wiki.orc.gmu.edu/mkdocs/Applying_for_an_Allocation/

### 7.3 Estimated Costs for This Project

| Run Type | Mesh Size | Cores | Walltime | SUs | Notes |
|----------|-----------|-------|----------|-----|-------|
| **Current RANS (laptop)** | 3.16M | 10 | 4 hrs | — | Free (WSL) |
| **RANS fine mesh** | 8M | 64 | 16 hrs | 1,024 | Mesh independence |
| **RANS extra-fine** | 20M | 128 | 40 hrs | 5,120 | Highest RANS fidelity |
| **Parametric RANS (10 angles)** | 8M | 64 | 160 hrs | 10,240 | Full wind rose |
| **LES single case** | 20M | 128 | 168 hrs | 21,504 | 7 days runtime |
| **LES parametric (3 angles)** | 20M | 128 | 504 hrs | 64,512 | Research-grade |

**Budget Recommendation:**
- **Phase 1 (validation):** 10,000 SUs (fine mesh + 5 wind angles)
- **Phase 2 (publication):** 50,000 SUs (LES + full parametric)

---

## 8. Best Practices & Tips

### 8.1 Storage Management
```bash
# Check your usage
du -sh /scratch/$USER/*
quota -s

# Compress old results
cd /scratch/$USER/openfoam_cases
tar -czf citicorp_3M_archive_$(date +%F).tar.gz citicorp_cfd/
# Then delete original: rm -rf citicorp_cfd/

# Scratch auto-purges after 30 days — transfer important results!
```

### 8.2 Checkpointing Long Jobs
```cpp
// In system/controlDict (for LES runs)
writeControl    adjustableRunTime;
writeInterval   10;          // Write every 10 seconds (adjust for storage)
purgeWrite      10;          // Keep only last 10 timesteps (save space)

// Restart if job times out:
startFrom       latestTime;  // Resume from last written time
```

### 8.3 Job Arrays for Parametric Studies
```bash
#!/bin/bash
#SBATCH --array=0-9             # 10 jobs: wind angles 0°, 10°, ..., 90°

# Calculate wind angle for this job
ANGLE=$(($SLURM_ARRAY_TASK_ID * 10))
echo "Running wind angle: ${ANGLE}°"

# Copy template case
cp -r /scratch/$USER/citicorp_template /scratch/$USER/citicorp_angle_${ANGLE}
cd /scratch/$USER/citicorp_angle_${ANGLE}

# Modify inlet direction (rotate U vector)
RAD=$(echo "scale=10; $ANGLE * 3.14159265 / 180" | bc)
UX=$(echo "scale=6; 44.7 * c($RAD)" | bc -l)
UY=$(echo "scale=6; 44.7 * s($RAD)" | bc -l)

sed -i "s/internalField   uniform (44.7 0 0)/internalField   uniform ($UX $UY 0)/" 0/U

# Run simulation
mpirun -np 64 simpleFoam -parallel > log.simpleFoam 2>&1
```

**Submit:**
```bash
sbatch job_parametric_array.slurm
# This creates 10 independent jobs, running in parallel if resources available
```

---

## 9. Quick Start Checklist

- [ ] **Account created:** https://wiki.orc.gmu.edu/mkdocs/Applying_for_an_Account/
- [ ] **SSH access working:** `ssh YOUR_GMU_ID@hopper.orc.gmu.edu`
- [ ] **Case transferred:** `rsync` or `scp` from WSL to `/scratch/$USER/`
- [ ] **OpenFOAM module loaded:** `module load openfoam/v2312`
- [ ] **Slurm script created:** Copy `job_rans_3M.slurm` above
- [ ] **Updated decomposeParDict:** `numberOfSubdomains` matches `#SBATCH --ntasks`
- [ ] **Job submitted:** `sbatch job_rans_3M.slurm`
- [ ] **Monitor:** `squeue -u $USER` and `tail -f log.simpleFoam`
- [ ] **Transfer results back:** `rsync` or `scp` to WSL/OneDrive

---

## 10. Example Session Transcript

```bash
# === On your laptop (WSL2) ===
cd ~/citicorp_cfd
git add -A && git commit -m "Ready for Hopper" && git push

# === SSH to Hopper ===
ssh YOUR_GMU_ID@hopper.orc.gmu.edu

# === On Hopper login node ===
cd /scratch/$USER
mkdir -p openfoam_cases
cd openfoam_cases

git clone https://github.com/YOUR_USERNAME/citicorp-crisis-simulation.git
cd citicorp-crisis-simulation/citicorp_cfd

# Load OpenFOAM
module load openfoam/v2312 impi/2021.5
source $FOAM_BASH

# Clean case
rm -rf processor* postProcessing log.*

# Update for 64 cores
sed -i 's/numberOfSubdomains.*/numberOfSubdomains 64;/' system/decomposeParDict

# Create job script (copy from Section 4.1 above)
nano job_rans_64cores.slurm
# (paste content, save with Ctrl+O, exit with Ctrl+X)

# Submit!
sbatch job_rans_64cores.slurm

# Check status
squeue -u $USER

# Monitor (wait until status = R for Running)
tail -f slurm_*.out

# When finished, check results
ls -lh 3000/
cat log.forces

# Transfer back
exit  # Logout from Hopper

# === Back on WSL2 ===
rsync -avz YOUR_GMU_ID@hopper:/scratch/$USER/openfoam_cases/citicorp-crisis-simulation/citicorp_cfd/3000/ \
    ~/citicorp_cfd_hopper_results/
```

---

## Additional Resources

- **Hopper Wiki:** https://wiki.orc.gmu.edu/mkdocs/
- **Slurm Cheat Sheet:** https://slurm.schedmd.com/pdfs/summary.pdf
- **OpenFOAM HPC Guide:** https://openfoamwiki.net/index.php/HowTo_running_OpenFOAM_in_parallel
- **GMU ORC Help:** orc@gmu.edu

---

**Next:** See [IMPROVEMENT_ROADMAP.md](IMPROVEMENT_ROADMAP.md) for technical details on improving the simulation fidelity.
