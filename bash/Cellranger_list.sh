#!/bin/bash -l

#SBATCH --job-name=cellranger # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 # CPU cores requested per task (change for threaded jobs)
#SBATCH --mem=128G # Memory needed per node (total)
#SBATCH --output=cellranger_%A_%a.out # File to which STDOUT will be written, including job ID

#---------------------Variables to be set-------------------------#
echo "Job ID : $JOB_ID"  ${SLURM_ARRAY_TASK_ID}
PROJECT_NAME="scRNAseq-Lung"
path=/athena/elementolab/scratch/yah2014/Projects/${PROJECT_NAME}
fastq_path=${path}/data/scRNA-seq/fastq
transcriptome="/athena/elementolab/scratch/yah2014/Indexed_genome/refdata-gex-GRCh38-2020-A"
fastq_dir_list="WC_30_T    WC_36_T"
sample_list=($fastq_dir_list)
sample_id=${sample_list[${SLURM_ARRAY_TASK_ID}]}
fastq_dir=Sample_$sample_id
echo $(ls -l $fastq_path/$fastq_dir)
localcores=16
localmem=128

cd $path/data/scRNA-seq/counts
# count
# Repeat this command per sample

cellranger_cmd="
cellranger count \
--id="${sample_id}" \
--sample="${sample_id}" \
--fastqs=$fastq_path/$fastq_dir  \
--transcriptome="${transcriptome}" \
--localcores="${localcores}" \
--localmem="${localmem}" \
--nosecondary \
--chemistry="auto"
"
echo -e "\n CMD: $cellranger_cmd \n"
$cellranger_cmd


