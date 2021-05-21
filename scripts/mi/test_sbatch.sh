#!/bin/bash
#SBATCH -J ecmb_test_sbatch
#SBATCH -o ./%N.%j.%a.out
#SBATCH -e ./%N.%j.%a.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=e.combrisson@gmail.com
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load all
module load anaconda/3
source activate ecmb_py38

python -O /home/combrisson.e/toolbox/seeg-ebrains/scripts/mi/compute_mi_nodes.py
