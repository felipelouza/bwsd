#! /bin/sh

#sources
K=(	10000	50000	100000	150000	200000	250000	300000
		10000000	50000000	100000000	150000000	200000000	250000000	300000000
		1000	5000	10000	50000	100000	150000	200000)

D=(	"all_ests.fasta" 
		"reads.fastq" 
		"citrus_ests.fasta")

##
host="jau"
#dir="/home/louza/dataset"
#workspace="/home/louza/Dropbox/workspace/"
dir="/mnt/data/bio/projects/gpt/"
workspace="/home/gpt/saca/"
exp="time"
result=1
##

#######################################################################

cd ${workspace}all-bwsd/

mkdir -p experiments/results 
make clean
make OUTPUT=${result}

for d in {0..2}
do

	test=${D[$d]}

	echo "${test}"

	path="experiments/results/${test}/"
	mkdir -p ${path} 

	for k in {0..6}
	do

		echo ${K[$k]}

		a=$(($d*7+$k))

		output=${path}${host}.${exp}.$k."txt"
		date >> ${output}

		echo $output
	
		for j in {1..3}	#MODE 1, 2 and 3
		do
	              
			echo " " >> ${output} 
			make run DIR=${dir} INPUT=${test} K=${K[$k]} MODE=$j >> ${output} 

		done
	
	done

done

#########################################################################

