#!/bin/sh
## Specify the name of your conda environment
SPK_ENV=spktest
## Specify the number of OpenMP threads
OMP_NUM_THREADS=8
## Specify the full path to conda
CONDA=conda

if ! command -v $CONDA > /dev/null 2>&1; then
    echo "Error: Conda command not found. Please install Conda or specify the full path in env.sh"
    echo "Visit https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html for installation instructions."
    exit 1
fi

eval "$($CONDA shell.bash hook)"
conda activate ${SPK_ENV}
export OMP_NUM_THREADS
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CONDA_PREFIX}/lib
