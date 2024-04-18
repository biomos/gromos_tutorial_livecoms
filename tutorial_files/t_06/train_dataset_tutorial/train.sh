#!/bin/bash

datapath='./meoh_trial.db' # your path to the ASE database
modelpath='./trial_model' # your path to the model directory

# model training
python spk_run.py train schnet custom $datapath $modelpath --property energy --derivative forces --rho property=0.01 derivative=0.99 --split 688 86 --batch_size 8 --n_epochs 2 --lr 0.0001 --lr_patience 10 --lr_min 1e-06 --cutoff 100.0 --num_gaussians 50 --features 32 --interactions 1

# model evaluation
python spk_run.py eval $datapath $modelpath --split test --batch_size 1
