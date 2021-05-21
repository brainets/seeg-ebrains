#!/bin/bash

# activate python
module purge
module load all
module load anaconda/3
source activate ecmb_py38

# run script
python -O /home/combrisson.e/toolbox/seeg-ebrains/scripts/mi/compute_mi_nodes.py
