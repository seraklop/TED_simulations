#!/bin/bash

nbTAXA=12
nbREPs=20
nbFACTORs=3
b=0

FOLDER=../a_simulateCharacters/Data_files
DIRNAME=${nbTAXA}taxa_allReps
mkdir $DIRNAME

for FAC in $( seq $nbFACTORs )
do	
	for REP in $( seq $nbREPs )
	do
		cp ${FOLDER}/Combined_b${b}_${nbTAXA}Taxa_${FAC}_${REP}.nex $DIRNAME
	done
done

cd ${DIRNAME}

cp ../c_Qsub_MB.sh .
cp ../d_input_${nbTAXA}Taxa.nex input.nex
sbatch c_Qsub_MB.sh

cd ..