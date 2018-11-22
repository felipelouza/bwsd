# vim: noai:ts=2:sw=2
#! /bin/sh

K=(	1500	3000	4500	6000	7500	9000	10500	12000	13500	150	17500	20000)

D=(	"ests.fasta" 
		"reads.fastq" 
		"uniprot.fasta"
		"wikipedia.txt")

#threads
T=(	1 2 4 8 16 32)

##
host="home"
dir="experiments/dataset/"
workspace="../"
exp="time"
result=1
##

#######################################################################

tar xvfz dataset.tar.gz -C ./

cd ${workspace} 
pwd

mkdir -p experiments/results 

#strings
for k in 0 #{0..9} 
do

echo ""
echo "####"

	#datasets
	for d in {0..3}
	do

		test=${D[$d]}
		echo "${test}"

		path="experiments/results/${test}/"
		mkdir -p ${path} 

		a=$k

		echo ${K[$a]}

		#n. of threads
		for t in 0 #{0..5}
		do

			echo ""
			echo "####"
			echo "Threads:" ${T[$t]}
			echo ""

			output=${path}${host}.${exp}.$a."t_"${T[$t]}."txt"
			date >> ${output}
			#echo $output
		
##
			for j in 3 1 2	
			do
				echo " " >> ${output} 
				rm -f ${dir}sdsl/*
				./bwsd ${dir}${test} ${K[$a]} -A $j -T ${T[$t]}  -v >> ${output} #-c
			done

##

			make WT=1 &>> out
			for j in 1 #MODE: wavelet-tree
			do
				echo " " >> ${output} 
				rm -f ${dir}sdsl/*
				./bwsd ${dir}${test} ${K[$a]} -A $j -T ${T[$t]}  -v >> ${output}
			done
##	
	
			make SD_VECTOR=0 &>> out
			for j in 1 #MODE: uncompressed bitvectors
			do
				echo " " >> ${output} 
				rm -f ${dir}sdsl/*
				./bwsd ${dir}${test} ${K[$a]} -A $j -T ${T[$t]}  -v >> ${output}
			done
	
##	
			for j in 4 #MODE: Alg. 2 in O(n+z) time
			do
				echo " " >> ${output} 
				rm -f ${dir}sdsl/*
				./bwsd ${dir}${test} ${K[$a]} -A $j -T ${T[$t]}  -v >> ${output}
			done
##	

		done

	done

done

#########################################################################

