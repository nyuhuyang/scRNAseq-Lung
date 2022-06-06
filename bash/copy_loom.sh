#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=velocyto_%A_%a.txt
#SBATCH --job-name=velocyto
#SBATCH --mem=2G  # memory requested, units available: K,M,G,T

echo $SLURM_ARRAY_TASK_ID
#---------------------Variables to be set-------------------------#
PROJECT_NAME="scRNAseq-AIM"
path=/athena/elementolab/scratch/yah2014/Projects/${PROJECT_NAME}/data
file_folder=$(ls ${path}/counts | tail -n +${SLURM_ARRAY_TASK_ID}| head -1) # Uses job array for each sample in the folder
#file="${file_folder}.bam" # add .bam
rmsk_gtf=/athena/elementolab/scratch/yah2014/Indexed_genome/hg38_rmsk.gtf
genes_gtf=/athena/elementolab/scratch/yah2014/Indexed_genome/refdata-gex-GRCh38-2020-A/genes/genes.gtf
echo "path="
echo "$path"
echo " "
echo $(ls -l $path/counts/$file_folder/)
echo $(ls -l $rmsk_gtf)
echo $(ls -l $genes_gtf)

echo "total size"
echo $(ls $path/counts/$file_folder)
echo " "
# ==== check if loom exist ======
echo $(ls -l $path/counts/$file_folder/velocyto/$file_folder.loom)
if test -f "$path/counts/$file_folder/velocyto/$file_folder.loom"
then
    echo "loom file exists"
    rsync --remove-source-files -av $path/counts/$file_folder/velocyto/$file_folder.loom $path/velocyto/$file_folder.loom
else
    echo "loom file does not exist."
fi
