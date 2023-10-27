#!/bin/bash
#SBATCH --job-name=dorado
#SBATCH --partition=owners,gpu
#SBATCH --gpus 1
#SBATCH --time=16:00:00
#SBATCH --mem=24G
#SBATCH -C GPU_GEN:VLT
#SBATCH -c 2
#SBATCH --chdir=/scratch/users/dubocd/altemose/slurm
#SBATCH --export=all
#SBATCH --requeue

source /home/users/dubocd/.bashrc
source activate general

export HDF5_PLUGIN_PATH=/scratch/users/dubocd/altemose/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin/


pod5=$1 # dir of pod5 files 
read_ids=$2


/oak/stanford/groups/altemose/dubocd/tools/dorado-0.3.4-linux-x64/bin/dorado basecaller \
--modified-bases 6mA 5mC \
-c 7500 \
-l $read_ids \
--emit-moves \
/oak/stanford/groups/altemose/dubocd/tools/dorado-0.3.4-linux-x64/dna_r10.4.1_e8.2_400bps_sup@v4.2.0 \
$pod5 \
> ${pod5}/dorado4_out.${read_ids}.bam 


exit

