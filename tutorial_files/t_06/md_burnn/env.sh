## Specify the name of your conda environment
SPK_ENV=spk
## Specify the number of OpenMP threads
OMP_NUM_THREADS=8

if ! command -v conda &> /dev/null; then
    echo "Error: Conda command not found. Please install Conda or ensure it is added to your PATH."
    echo "Visit https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html for installation instructions."
    exit 1
fi

eval "$(conda shell.bash activate)"
conda activate ${SPK_ENV}
export OMP_NUM_THREADS
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CONDA_PREFIX}/envs/${SPK_ENV}/lib
