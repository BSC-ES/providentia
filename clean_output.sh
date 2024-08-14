#!/bin/bash
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

rm -rf figures/*
rm -rf plots/*
rm -rf reports/*
rm -rf saved_data/*
rm -rf providentia/__pycache__
rm -rf notebooks/.ipynb_checkpoints
rm -f interactive.out

# Interpolation
rm -rf providentia/interpolation/arguments/*
rm -rf providentia/interpolation/interpolation_logs/*
rm -rf providentia/interpolation/management_logs/*
rm -rf providentia/interpolation/submission_logs/*
rm -rf providentia/interpolation/submit/*
rm -rf providentia/interpolation/__pycache__ 
rm -rf providentia/interpolation/slurm*