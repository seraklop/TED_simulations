#!/bin/sh

module add vital-it
module add R/latest

nbTAXA=12
nbREPs=20
nbFACTORs=3
b=0

outFILE="b${b}_Results_${nbTAXA}taxa.summary"
echo "Rep Fac lnLharmMean ASDSF RootMean RootMedian RootLower RootUpper RootPSRF m1Mean m1Median m1Lower m1Upper m1PSRF m2Mean m2Median m2Lower m2Upper m2PSRF CRMean CRMedian CRLower CRUpper CRPSRF propCorr_wFos propCorr_woFos" > ${outFILE}

cd ${nbTAXA}taxa_allReps
cp ../f_calculate_topoDist.r calculate_topoDist.r

for REP in $( seq $nbREPs )
do
	echo "analyzing data in ${nbTAXA}taxa_rep_${REP} "	
	sed -n "${REP}p" ../Trees_${nbTAXA}Taxa_20.tre > reftree.nwk
	
	for DATA in *_${REP}.nex.lstat
	do
		
		fileSTEM=$( echo ${DATA} | cut -d "." -f1 ).nex
		FAC=$( echo ${fileSTEM} | cut -d "_" -f4 )
		lnL=$( grep all ${DATA} | cut -f3 )
		ASDSF=$( awk '{print $NF}' ${fileSTEM}.mcmc | tail -n 1 )
		RMean=$( grep 'age{all}\[0\]' ${fileSTEM}.vstat | cut -f2 )
		RMedian=$( grep 'age{all}\[0\]' ${fileSTEM}.vstat | cut -f6 )
		RLower=$( grep 'age{all}\[0\]' ${fileSTEM}.vstat | cut -f4 )
		RUpper=$( grep 'age{all}\[0\]' ${fileSTEM}.vstat | cut -f5 )
		RPSRF=$( grep 'age{all}\[0\]' ${fileSTEM}.vstat | cut -f7 )
		m1Mean=$( grep 'm{1}' ${fileSTEM}.pstat | cut -f2 )
		m1Median=$( grep 'm{1}' ${fileSTEM}.pstat | cut -f6 )
		m1Lower=$( grep 'm{1}' ${fileSTEM}.pstat | cut -f4 )
		m1Upper=$( grep 'm{1}' ${fileSTEM}.pstat | cut -f5 )
		m1PSRF=$( grep 'm{1}' ${fileSTEM}.pstat | cut -f9 )
		m2Mean=$( grep 'm{2}' ${fileSTEM}.pstat | cut -f2 )
		m2Median=$( grep 'm{2}' ${fileSTEM}.pstat | cut -f6 )
		m2Lower=$( grep 'm{2}' ${fileSTEM}.pstat | cut -f4 )
		m2Upper=$( grep 'm{2}' ${fileSTEM}.pstat | cut -f5 )
		m2PSRF=$( grep 'm{2}' ${fileSTEM}.pstat | cut -f9 )
		CRMean=$( grep 'lockrate' ${fileSTEM}.pstat | cut -f2 )
		CRMedian=$( grep 'lockrate' ${fileSTEM}.pstat | cut -f6 )
		CRLower=$( grep 'lockrate' ${fileSTEM}.pstat | cut -f4 )
		CRUpper=$( grep 'lockrate' ${fileSTEM}.pstat | cut -f5 )
		CRPSRF=$( grep 'lockrate' ${fileSTEM}.pstat | cut -f9 )
		
		cp ${fileSTEM}.con.tre con.tre	
		cp ${fileSTEM}.withoutFossils.con.tre con_woF.tre
		R CMD BATCH --slave calculate_topoDist.r
		
		propCorr_wFOS=$( grep 'withFossils' nodes_correct.txt | cut -f2 )
		propCorr_woFOS=$( grep 'woFossils' nodes_correct.txt | cut -f2 )			
		rm con.tre
		rm con_woF.tre
		rm nodes_correct.txt
		
		echo ${REP} ${FAC} $lnL $ASDSF $RMean $RMedian $RLower $RUpper $RPSRF $m1Mean $m1Median $m1Lower $m1Upper $m1PSRF $m2Mean $m2Median $m2Lower $m2Upper $m2PSRF $CRMean $CRMedian $CRLower $CRUpper $CRPSRF $propCorr_wFOS $propCorr_woFOS >> ../${outFILE} 
	done
	rm reftree.nwk		
	
done

