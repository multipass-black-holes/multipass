#!/bin/bash
#SBATCH --partition=shared
#SBATCH --time=2-00:00:00

#SBATCH --ntasks=30
#SBATCH --error=%x/long-%A.err
#SBATCH --output=%x/long-%A.out
#SBATCH --mem=50G

module load devel/CMake/3.24.3-GCCcore-12.2.0
module load compiler/GCC/12.2.0
module load lang/Python/3.10.8-GCCcore-12.2.0
module load toolchain/foss/2022b

#srun --mpi=pmix -n ${SLURM_NTASKS} ./main $SLURM_JOB_NAME/long.card
mpirun -n ${SLURM_NTASKS} ./main $SLURM_JOB_NAME/long.card
