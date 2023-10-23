#!/bin/bash
#SBATCH --job-name=mm2
#SBATCH --partition=owners
#SBATCH --time=16:00:00
#SBATCH --mem=100G
#SBATCH -c 20
#SBATCH --chdir=/scratch/users/dubocd/altemose/slurm
#SBATCH --export=all
#SBATCH --requeue

source /home/users/dubocd/.bashrc
source activate general

reads=$1
ref=$2

ref_basename=$(basename $ref)

minimap2 -ax map-ont --MD -Y --eqx -t 20 $ref $reads | samtools sort -O BAM -@ 20 > ${reads}.aligned_${ref_basename}.bam 

exit
