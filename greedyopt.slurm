#!/usr/bin/env bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -t 10:00:00
#SBATCH -J greedyoptjob
#SBATCH -o greeydoptjob.out
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=brush@princeton.edu

echo '-------------------------------'
${SLURM_SUBMIT_DIR}
echo ${SLURM_SUBMIT_DIR}
echo Running on host $(hostname)
echo Time is $(date)
echo SLURM_NODES are $(echo ${SLURM_NODELIST})
echo '-------------------------------'
echo -e '\n\n'

EXEC=/usr/licensed/bin/matlab
OPTS=' -singleCompThread -nodisplay -nosplash'

export PROCS=${SLURM_CPUS_ON_NODE}

${EXEC}${OPTS} < greedy_opt_v2.m