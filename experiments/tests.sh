#! /bin/sh

#sources
K=(	1500	3000	4500	6000	7500	9000	10500	12000	13500	15000 	17500	20000)
#K=(	1500	3000	4500	6000	7500	9000	10500
#	2500	5000	7500	10000	12500	15000	17500
#	1000	2000	3000	4000	5000	6000	7000	
#	150	300	450	600	750	900	1050)

D=(	"all_ests.fasta" 
	"reads.fastq" 
	"uniprot_trembl-201509.fasta"
	"wikipedia-SURF.txt")

##
host="jau"
#dir="/home/louza/dataset"
#workspace="/home/louza/Dropbox/workspace/"
dir="/mnt/data/bio/projects/gpt/text/"
workspace="/home/gpt/saca/"
exp="time"
result=1
##

#######################################################################

#cd ${workspace}all-bwsd/
cd ..

mkdir -p experiments/results 

for k in {0..9} #7 8 9 0 1 2 3 4 5 6 10 11 #{0..6} 
do

echo "####"

	for d in {0..3}
	do

		test=${D[$d]}

		echo "${test}"

		path="experiments/results/${test}/"
		mkdir -p ${path} 


		#a=$(($d*7+$k))
		a=$k

		echo ${K[$a]}

		output=${path}${host}.${exp}.$a."txt"
		date >> ${output}

		echo $output
	

###
#		rm -f ${dir}sdsl/*
#		make OUTPUT=${result} WT=1 &>> out
#
#		for j in 3 1 #{1..3}	#MODE: wavelet-tree
#		do
#			echo " " >> ${output} 
#			make run DIR=${dir} INPUT=${test} K=${K[$a]} MODE=$j >> ${output} 
#
#		done
###	
#
#		rm -f ${dir}sdsl/*
#		make OUTPUT=${result} SD_VECTOR=0 &>> out
#
#		for j in 1 #MODE: bitmap
#		do
#			echo " " >> ${output} 
#			make run DIR=${dir} INPUT=${test} K=${K[$a]} MODE=$j >> ${output} 
#
#		done
###	
#
#		rm -f ${dir}sdsl/*
#		make OUTPUT=${result} &>> out
#
#		for j in 1 #MODE: bitmap_sd
#		do
#			echo " " >> ${output} 
#			make run DIR=${dir} INPUT=${test} K=${K[$a]} MODE=$j >> ${output} 
#
#		done
##	

		rm -f ${dir}sdsl/*
		make OUTPUT=${result} &>> out

		for j in 2 #MODE: Gog's
		do
			echo " " >> ${output} 
			make run DIR=${dir} INPUT=${test} K=${K[$a]} MODE=$j >> ${output} 

		done
##	


	done

done

#########################################################################

