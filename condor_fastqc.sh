#! /bin/bash

# This script assums a project directory structure
# 	project_name
#		|__data
#			|__file1.fq
#			|__file2.fq	
#			|....
# 	 	|__results
#  	|__scripts
#			|__<scripts for further analysis>
#
# and will result in the directory structure
# 	project_name
#			|__file1.fq
#			|__file2.fq	
#			|....	
# 	 	|__results
#			|__QC_<timestamp>
#				|__condor_qc_<timestamp>
#				|__logs
#					|__*.err
#					|__*.out
#					|__*.log
#					|__...(etc. for other jobs)
#				|__<fastq>.html
#				|__<fastq>.zip
#				|__<fastq>
#				|__... (and so on for the other fasta files)
#  	|__scripts
#			|__<scripts for further analysis>
#
# This script can be run from inside the "project_name" directory, or
# pass in the "project_name" directory as an argument. When called
# this script generates a condor submit file called "condor_qc_<timestamp>"
# When called, it will create a condor submit file
# named condor_submit_file_<timestamp> which is placed
# in 'scripts'. This script then submits that file
# to the condor queue. The jobs will run fastqc on each fasta
# file in the data directory and put the results in a directory
# 'results/QC_<timestamp>' .
#
# The only required argument is the file extension for the input files.
# The project directory can be passed in as an option argument if notification
# run in the project directory.

file_extension=$1
project_directory=$2

project_directory=${project_directory:=$(pwd -P)}
echo "project directory: $project_directory"

data_directory=$project_directory/data
echo "data directory: $data_directory"

results_directory=$project_directory/results
scripts_directory=$project_directory/scripts

timestamp=$(date | sed 's/ /_/g' | sed 's/\:/-/g')
output_directory=$results_directory/QC_$timestamp
#output_directory=$results_directory/QC
echo "output directory: $output_directory"

echo `mkdir $output_directory`
echo `mkdir $output_directory/logs`

common_dir=$project_directory
request_cpus=1
request_memory=56GB
email=balter@ohsu.edu
getenv="True"
notification=Error
max_execution_time=10000 

submit_file=condor_qc_$timestamp
#submit_file=condor_qc 

printf "
executable				=/mnt/lustre1/CompBio/bin/fastqc
output						=$output_directory/logs/\$(Cluster).\$(Process).out
error							=$output_directory/logs/\$(Cluster).\$(Process).err
log							=$output_directory/logs/\$(Cluster).\$(Process).log
request_cpus				=$request_cpus
request_memory		=$request_memory
notify_user				=$email
notification				=$notification
getenv						=$getenv
+MaxExecutionTime	=$max_execution_time

" \
> ${submit_file}

#for fasta_file in $(ls $data_directory *.fq.gz)
#for fasta_file in $(find data/fasta_files/ -type f -name "*.fq.gz" -exec basename {} \;)
for fasta_file in $data_directory/*.$file_extension
do
	echo $(basename $fasta_file)
	echo "arguments=$fasta_file --outdir=$output_directory --extract"
	printf "
arguments=$fasta_file --outdir=$output_directory --extract
queue 1
" \
>>${submit_file}
done

echo `mv $submit_file $output_directory`
echo "submit file: $submit_file"

echo `condor_submit $output_directory/$submit_file`

#cat submit_file
