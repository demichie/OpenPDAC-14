#!/bin/bash -l
#
#SBATCH --cluster=ub-hpc
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=100000
#SBATCH --job-name="ValenSingleP"
#SBATCH --output=log2.foamRun
#SBATCH --mail-type=all

module load gcc/11.2.0
module load openmpi/4.1.1
module load trilinos/13.4.1
module load scotch/7.0.2

srun -n 10 foamRun -parallel

reconstructPar -noLagrangian > log.reconstruct
rm -rf processor*
