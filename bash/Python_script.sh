#!/bin/bash -l

#SBATCH --partition=panda   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=scRNAseq-velocyto.txt
#SBATCH --job-name=scRNAseq-velocyto
#SBATCH --mem=80G   # memory requested, units available: K,M,G,T


# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID"  $SLURM_ARRAY_TASK_ID"
echo "=========================================================="

file=/athena/elementolab/scratch/yah2014/Projects/scRNAseq-Lung/bash/velocyto_single.py
echo $(ls -l $file)
python $file $SLURM_ARRAY_TASK_ID
