#!/bin/bash
#PBS -q standard
#PBS -l select=1:ncpus=1:mem=6gb
#PBS -N bATAC_batch
#PBS -W group_list=cusanovichlab
#PBS -l place=pack:shared
#PBS -l walltime=2:00:00
### cput ("total cputime") = walltime * ncpus
##PBS -l cput=16:00:00 

# need to change fastq dir in line 14 and script name in line 21

date
for file in /rsgrps/cusanovichlab/projects/bulk_atac/various_060719/fastq/merged.fastq/*R1.fastq.gz
do
	base=$(basename ${file} _R1.fastq.gz)
	fastq1=${base}_R1.fastq.gz
	fastq2=${base}_R2.fastq.gz
	#Use double quotes to make the shell expand variables while preserving whitespace
	#number before s is the line number to search; number at the end of sed expression is the number of occurrence to replace
	sed "5 s/bATAC/$base/" /home/u23/haozhang1/bATAC_pipeline_single_job_mm.sh > /home/u23/haozhang1/bATAC_script_batch/bATAC_pipeline_${base}.sh
	sed -i "26 s/fastq1/$fastq1/2" /home/u23/haozhang1/bATAC_script_batch/bATAC_pipeline_${base}.sh
	sed -i "27 s/fastq2/$fastq2/2" /home/u23/haozhang1/bATAC_script_batch/bATAC_pipeline_${base}.sh
	sed -i "28 s/base/$base/2" /home/u23/haozhang1/bATAC_script_batch/bATAC_pipeline_${base}.sh
	# /^word$/ to concatenate the next line only when your current line is word
	#N; concatenate two lines
	#unfortunately, not work on UA HPC
	#sed -i "/^fastq1=fastq1$/ {N; s/fastq1=fastq1\nfastq2=fastq2/fastq1=$fastq1\nfastq2=$fastq2/}" /home/u23/haozhang1/bATAC_script_batch/bATAC_pipeline_${base}.sh
	qsub /home/u23/haozhang1/bATAC_script_batch/bATAC_pipeline_${base}.sh
done
date

