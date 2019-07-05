#!/bin/bash

#This script takes PHYLIP input files, 1 for molecular partition (a single alignment in the file)
#and X for morphology partition (with as many datasets in it as you like, but all the same size), 
#and combines them into X NEXUS matrices for MrBayes

#Change the settings in these first lines
REPS=20
NbTAXA=12

#--------------------------------
#main part
#--------------------------------

for (( i=1; i<=$REPS; i++))
do
	INFILE_MOLEC=Data_files/Seqs_${NbTAXA}Taxa_${i}.phy  
	INFILE_MORPH=Data_files/MorphoData_${NbTAXA}Taxa_${i}.phy  

	#file conversion - as unix an windows/mac have different linebreak characters
	tr -d '\r' < $INFILE_MOLEC > temp
	mv temp $INFILE_MOLEC
	tr -d '\r' < $INFILE_MORPH > temp
	mv temp $INFILE_MORPH

	#read the numbers of taxa and sites from file
	NBTAX=$(head -n 1 $INFILE_MOLEC | cut -d " " -f 1)
	NBTAX2=$(head -n 1 $INFILE_MORPH | cut -d " " -f 1)

	if [ $NBTAX -ne $NBTAX2 ]; then
   		echo "Error: different nb of taxa in the molecular and morpho file. exiting..."
   		exit 1
	fi

	NBSITES=$(head -n 1 $INFILE_MOLEC | cut -d " " -f 2)
	NBCHARS=$(head -n 1 $INFILE_MORPH | cut -d " " -f 2)
	let NBTOT="$NBSITES + $NBCHARS"
	let DNASTART="$NBCHARS + 1"

	#calculate the numbers of datasets based on the nb of lines in the morpho file and the nb of taxa
	NBLINES=$(wc -l $INFILE_MORPH | cut -d " " -f 1)
	let NBDATASETS="($NBLINES)/($NBTAX+2)"

	if [ $NBDATASETS -eq 0 ]; then
		NBDATASETS=1
	fi

	for (( c=1; c<=$NBDATASETS; c++))
	do
		#write header of NEXUS file including matrix size
		X=$(printf %04d ${c%})
		OUTFILENAME=Data_files/Combined_${NbTAXA}Taxa_${i}_$X.nex
		echo "#NEXUS" > $OUTFILENAME
		echo "begin data;" >> $OUTFILENAME
		echo dimensions ntax=$NBTAX nchar=$NBTOT\; >> $OUTFILENAME	
	 	echo format datatype=mixed\(standard:1-$NBCHARS,DNA:$DNASTART-$NBTOT\) interleave=yes gap=- missing=?\; >> $OUTFILENAME
	 	echo matrix >> $OUTFILENAME
	 	
		#now get the current line numbers where the morpho dataset is in the morpho file. add
		let MORPHSTART="($c-1)*($NBTAX+2)+2"
		let MORPHEND="($c)*($NBTAX+2)"
		
		for (( l=$MORPHSTART; l<=$MORPHEND; l++))
		do
			sed -n ${l}p $INFILE_MORPH >> $OUTFILENAME
		done
	
		#and now we also add the molecular data below, and end the NEXUS block
		let MOLEND="$NBTAX +1"
		for (( m=2; m<=$MOLEND; m++))
		do
			sed -n ${m}p $INFILE_MOLEC >> $OUTFILENAME
		done
		echo \; >> $OUTFILENAME
		echo end\; >> $OUTFILENAME
	done
done

#-------------
