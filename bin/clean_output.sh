#!/bin/bash
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

rm -rf ./figures/*
rm -rf ./plots/*
rm -rf ./reports/*
rm -rf ./saved_data/*
rm -rf ./providentia/__pycache__
rm -rf ./notebooks/.ipynb_checkpoints
rm -f ./notebook.out

# Interpolation
rm -rf ./logs/interpolation/arguments/*
rm -rf ./logs/interpolation/interpolation_logs/*
rm -rf ./logs/interpolation/management_logs/*
rm -rf ./logs/interpolation/submission_logs/*
rm -rf ./logs/interpolation/greasy_logs/*
rm -rf ./logs/interpolation/__pycache__ 
rm -rf ./logs/interpolation/slurm*