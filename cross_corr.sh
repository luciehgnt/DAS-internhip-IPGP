#!/bin/bash
#SBATCH --job-name CROSS_CORR
#SBATCH --output CROSS_CORR_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --partition ncpu,ncpum

echo "Lanzado en: $SLURM_NODELIST"
echo "SLURM_NTASKS_PER_NODE=$SLURM_NTASKS_PER_NODE"
echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"
echo "SLURM_NNODES=$SLURM_NNODES"
echo "SLURM_CPUS_ON_NODE=$SLURM_CPUS_ON_NODE"

cd $SLURM_SUBMIT_DIR

# module purge
# module load openmpi/gcc-10.3.0/4.1.1
# module load gcc/10.3.0
# module load cuda/11.2
#module load intel/21U2/suite
#module load anaconda3

export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
# export OMP_SCHEDULE="static,1"
# export OMP_DYNAMIC="false"
# export OMP_PROC_BIND="true"
# export OMP_NESTED="false"
# export OMP_MAX_ACTIVE_LEVELS="2"

# export OMP_WAIT_POLICY="active"
# export OMP_STACKSIZE="32k"
# export OMP_THREAD_LIMIT="128"

echo "$OMP_NUM_THREADS"
eval "$(conda shell.bash hook)"
conda activate das
echo "Go !!"
python3 cross_corr.py