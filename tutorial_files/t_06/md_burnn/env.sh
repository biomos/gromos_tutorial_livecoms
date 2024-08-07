## Specify the name of your conda environment
SPK_ENV=spk
## Specify the number of OpenMP threads
OMP_NUM_THREADS=8


. ${CONDA_PREFIX}/etc/profile.d/conda.sh
conda activate ${SPK_ENV}
export OMP_NUM_THREADS
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CONDA_PREFIX}/envs/${SPK_ENV}/lib
