#!/bin/bash
#module load nvhpc_sdk/20.11

#SBATCH --partition=v100
#SBATCH --gres=gpu:2


#SBATCH --job-name=osctest
#SBATCH --output=osc.txt

#SBATCH --partition=short
#SBATCH --cpus-per-task=20
#SBATCH --ntasks-per-core=1

#nsys nvprof -t  --print-gpu-trace ./$program
#nsys profile --stats=true ./$program --id $id --conf $conf
#srun -p $target_machine --gres=gpu:$ngpu ./$program --id $id --conf $conf
