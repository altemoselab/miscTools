#!/bin/bash
#SBATCH --job-name=guppy
#SBATCH --partition=owners,gpu,sfgf
#SBATCH --gpus 1
#SBATCH --time=16:00:00
#SBATCH --mem=24G
#SBATCH -c 4
#SBATCH --chdir=/scratch/users/dubocd/altemose/slurm
#SBATCH --export=all
#SBATCH --requeue


source /home/users/dubocd/.bashrc
source activate general

ml gcc
ml cuda
ml cudnn 
ml nvidia-smi 

# full paths required for all args

fast5_dir=$1
output_dir=$2
model=$3 # .cfg file for model, needs to be in the guppy/data directory directory

/oak/stanford/groups/altemose/dubocd/rotation/guppy/ont-guppy/bin/guppy_basecaller \
-i $fast5_dir\
-s $output_dir \
--bam_out --bam_methylation_threshold 0 \
-c $model \
--device cuda:all 

exit

