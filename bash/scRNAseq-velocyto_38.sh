#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --output=velocyto_%A_%a.txt
#SBATCH --job-name=velocyto
#SBATCH --mem=64G  # memory requested, units available: K,M,G,T

conda activate r4.0.3
spack load -r samtools@1.8
samtools --version
echo $SLURM_ARRAY_TASK_ID
#---------------------Variables to be set-------------------------#
PROJECT_NAME="scRNAseq-Lung"
path=/athena/elementolab/scratch/yah2014/Projects/${PROJECT_NAME}/data/scRNA-seq
file_folder=$(ls ${path}/counts | tail -n +${SLURM_ARRAY_TASK_ID}| head -1) # Uses job array for each sample in the folder
#file="${file_folder}.bam" # add .bam
rmsk_gtf=/athena/elementolab/scratch/yah2014/Indexed_genome/hg38_rmsk.gtf
genes_gtf=/athena/elementolab/scratch/yah2014/Indexed_genome/refdata-gex-GRCh38-2020-A/genes/genes.gtf
echo "path="
echo "$path"
echo " "
echo $(ls $path/counts/$file_folder/)
echo $(ls -l $rmsk_gtf)
echo $(ls -l $genes_gtf)

# ==== check if loom exist ======
if test -f "$path/velocyto/$file_folder.loom"
then
    return 0;
    echo "loom file exist"
else
    echo "loom file does not exist. Continue velocyto run10x."
fi
#----------------sort BAM File-------------------
echo $(ls -l $path/counts/$file_folder/outs/possorted_genome_bam.bam)
if test -f "$path/counts/$file_folder/outs/cellsorted_possorted_genome_bam.bam"
then
    echo "sorted bam exists"
else
    echo "sorted bam doesn't exists"
    echo "to sort by cellID"
    samtools sort -t CB -O BAM -o $path/counts/$file_folder/outs/cellsorted_possorted_genome_bam.bam $path/counts/$file_folder/outs/possorted_genome_bam.bam
fi

#-----------velocyto Command--------------------------------#
echo "Processing velocyto run10x"
echo " "
echo "-------------------------------- "
echo "Processing $file_folder"
cd $path/velocyto
velocyto run10x -m $rmsk_gtf $path/counts/$file_folder $genes_gtf --samtools-memory 250000

# ==== check if loom exist ======
echo "velocyto output files:"
echo " "
echo $(ls -l $path/counts/$file_folder/velocyto/$file_folder.loom)
if test -f "$path/counts/$file_folder/velocyto/$file_folder.loom"
then
    echo "loom file exists"
    rsync --remove-source-files -av $path/counts/$file_folder/velocyto/$file_folder.loom $path/velocyto/$file_folder.loom
    echo " "
    echo "velocyto run10x Complished"
else
    echo "loom file does not exist."
    echo " "
    echo "velocyto run10x Failed."
fi

if test -f "$path/counts/$file_folder/outs/cellsorted_possorted_genome_bam.bam"
then
    echo "sorted bam exists"
    rm $path/counts/$file_folder/outs/cellsorted_possorted_genome_bam.bam*
fi
