#!/bin/bash
#SBATCH --job-name=dorado_v5.3
#SBATCH --partition=owners,gpu
#SBATCH --gpus 1
#SBATCH --time=16:00:00
#SBATCH --mem=48G
#SBATCH -C GPU_GEN:VLT
#SBATCH -c 8
#SBATCH --chdir=/scratch/users/dubocd/altemose/slurm
#SBATCH --export=all
#SBATCH --requeue

source /home/users/dubocd/.bashrc
source activate general

pod5=$1 # dir of pod5 files 
read_ids=$2
reference=$3

base_read_id_file=$(basename $read_ids)

/oak/stanford/groups/altemose/tools/dorado-0.5.3-linux-x64/bin/dorado basecaller \
--modified-bases 5mC_5hmC \
-l $read_ids \
--reference $reference \
/oak/stanford/groups/altemose/tools/dorado-0.5.1-linux-x64/models/dna_r10.4.1_e8.2_400bps_sup@v4.3.0 \
$pod5 > ${read_ids}.dorado_v5.3_mC.bam


exit

