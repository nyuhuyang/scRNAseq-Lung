#!/bin/bash -l

#SBATCH --partition=panda   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=STAR_%A_%a.txt
#SBATCH --job-name=STAR
#SBATCH -c 32 		# number of cores requested -- this needs to be greater than or equal to the number of cores you plan to use to run your job
#SBATCH --mem=64G   # memory requested, units available: K,M,G,T

echo "Job ID : $JOB_ID"  ${SLURM_ARRAY_TASK_ID}

#---------------------Variables to be set-------------------------#
conda activate rnaseq
PROJECT_NAME="scRNAseq-Lung"
path=/athena/elementolab/scratch/yah2014/Projects/${PROJECT_NAME}
STARREF=/athena/elementolab/scratch/yah2014/Indexed_genome/hg38_Genomedir/
gtf="/athena/elementolab/scratch/yah2014/Indexed_genome/reference_sources/gencode.v42.primary_assembly.annotation.gtf"
fastq_path=${path}/data/RNA-seq/fastq
fastq_dir_list="SR_EX42_CTR_S1    SR_EX42_IFN_S2    RS58_3_S3    RS58_5_S4"
dir_list=($fastq_dir_list)
fastq_dir=${dir_list[${SLURM_ARRAY_TASK_ID}]}
echo $(ls -l $fastq_path/$fastq_dir)

Sample=$fastq_dir
echo "path= $path"
echo "STARREF= $(ls -l $STARREF)"
echo "GTF= $(ls -l $gtf)"
echo "Sample= $Sample"

#mkdir ${path}/data/RNA-seq/BAMS
#mkdir ${path}/data/RNA-seq/FPKM
#mkdir ${path}/data/RNA-seq/Counts


#--------Formatting multiple lane Data for STAR-------------------#
F1=$(ls $fastq_path/$fastq_dir/*R1*)
F2=$(ls $fastq_path/$fastq_dir/*R2*)
F1=$(echo $F1|tr ' ' ',') #Turn space into ,
F2=$(echo $F2|tr ' ' ',') #Turn space into ,
echo "F1= " ${F1}
echo "F2= " ${F2}

#-----------STAR Alignment Command--------------------------------#
echo "Processing Sample"
echo " "
mkdir ${TMPDIR}/${Sample}
cd ${TMPDIR}/${Sample}
echo "-------------------------------- "
echo 'Processing' ${Sample}
echo "Aligning FASTQ"
STAR --genomeDir ${STARREF} \
--runMode alignReads \
--readFilesIn ${F1} ${F2} \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${Sample}_  \
--runThreadN 12 \
--readFilesCommand zcat \
--outFilterIntronMotifs RemoveNoncanonical
echo "STAR alignment Complished"
echo "STAR output files:"
echo $(ls -l )
echo " "
#----------------Sorting SAM File by name-------------------
echo "Samtools start"
#samtools view -bS *Aligned.out.sam > ${Sample}_Aligned.out.bam
#samtools sort -n -O BAM ${Sample}_Aligned.out.bam -o ${Sample}_name_sorted.bam #no need for new HTSeq
#samtools sort -O BAM ${Sample}_Aligned.out.bam -o ${Sample}_sorted.bam # for Cufflinks

#----------------Samtools index ------------------------------

#samtools index ${Sample}_sorted.bam
samtools index ${Sample}_Aligned.sortedByCoord.out.bam

rsync -r -v --exclude="*.sam" $TMPDIR/${Sample} ${path}/data/RNA-seq/BAMS/

echo $(ls -l )


#-------------------CuffLinks Command---------------------------
echo " "
echo "Cufflinks Start"
cufflinks -q --max-bundle-frags 100000000  \
-p 2 \
-G $gtf \
-o $TMPDIR/${Sample}_CuffLinks  $TMPDIR/${Sample}/${Sample}_Aligned.sortedByCoord.out.bam
echo "Cufflinks End"
rsync -r -v $TMPDIR/${Sample}_CuffLinks ${path}/data/RNA-seq/FPKM
#---------------------------------------------------------------

#--------------Running HTSeq Count-----------------------------
# a Running htseq-count
echo " "
echo "HTSeq Count Start"
htseq-count -s no -f bam $TMPDIR/${Sample}/${Sample}_Aligned.sortedByCoord.out.bam $gtf \
> $TMPDIR/${Sample}.bam.count

rsync -r -v $TMPDIR/${Sample}.bam.count  ${path}/data/RNA-seq/Counts
echo "HTSeq Count End"
