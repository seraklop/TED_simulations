#!/bin/sh

#SBATCH --mail-user=seraina.klopfstein@iee.unibe.ch
#SBATCH --mail-type=fail,end
#SBATCH --job-name="6.6b0_MCl"
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=all
#SBATCH --output=outfile
#SBATCH --error=errfile

module add vital-it/7
module add Phylogeny/mrbayes/3.2.6

for DATA in Combined_*.nex
do
	cp $DATA data.nex
	tr -d '\r' < input.nex > temp
	mv temp input.nex
	
	mb input.nex
	
	rm data.nex.run2.t
	rm data.nex.run2.p
	rm data.nex.ck*
	rm data.*.trprobs
	for FILE in data.nex.*
	do
		FILEext=$( echo $FILE | sed 's:data.nex.::' )
		mv "$FILE" "${DATA}.${FILEext}"
	done
done

