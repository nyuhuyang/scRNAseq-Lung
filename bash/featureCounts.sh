#!/bin/bash -l

#SBATCH --partition=panda   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=STAR_%A_%a.txt
#SBATCH --job-name=STAR
#SBATCH -c 32         # number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --mem=64G   # memory requested, units available: K,M,G,T

echo "Job ID : $JOB_ID"  ${SLURM_ARRAY_TASK_ID}

#---------------------Variables to be set-------------------------#
conda activate rnaseq
PROJECT_NAME="scRNAseq-Lung"
path=/athena/elementolab/scratch/yah2014/Projects/${PROJECT_NAME}
gtf=/athena/elementolab/scratch/yah2014/Indexed_genome/reference_sources/Homo_sapiens.GRCh38.108.gtf
Sample="SR_EX42_CTR_S1.SR_EX42_IFN_S2.RS58_3_S3.RS58_5_S4"
Sample_list=("SR_EX42_CTR_S1" "SR_EX42_IFN_S2"  "RS58_3_S3"    "RS58_5_S4")

#--------------Running featureCounts Count-----------------------------
# a Running featureCounts Count
echo " "
echo "featureCounts Count Start"
featureCounts -a $gtf -o ${path}/data/RNA-seq/Counts/${Sample}.bam.count -T 24 ${path}/data/RNA-seq/BAMS/*_Aligned.sortedByCoord.out.bam

picard CollectInsertSizeMetrics I=/athena/elementolab/scratch/yah2014/Projects/scRNAseq-Lung/data/RNA-seq/BAMS/RS58_3_S3_Aligned.sortedByCoord.out.bam O=insert_size_metrics.txt H=insert_size_histogram.pdf M=0.5
for sample in ${Sample_list[@]}; do picard CollectInsertSizeMetrics I=${path}/data/RNA-seq/BAMS/${sample}_Aligned.sortedByCoord.out.bam O=${path}/data/RNA-seq/Counts/${sample}_insert_size_metrics.txt H=${path}/data/RNA-seq/Counts/${sample}_insert_size_histogram.pdf M=0.5; done
