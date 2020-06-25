#!/bin/bash

## MAKE SURE YOU ARE ALREADY IN THE QIIME1 ENVIRONMENT
# Workflow to run all of the SALBIN at once

biomBALTIC=./INPUT_DATA/Baltic_16S/OTU_Table_final.txt
MFBALTIC=./INPUT_DATA/Baltic_16S/metadata_table.tsv
outputBALTIC=01_16sBaltic_QTAG

biomFRASER16=./INPUT_DATA/Fraser_16S/OTU_Table_final.txt
MFFRASER16=./INPUT_DATA/Fraser_16S/MF_16sFraser_noConCOL.txt
outputFRASER16=01_16sFraser_QTAG

biomFRASER18=./INPUT_DATA/Fraser_18S/OTU_Table_final.txt
MFFRASER18=./INPUT_DATA/Fraser_18S/MF_18sFraser_noConCOL.txt
outputFRASER18=01_18sFraser_QTAG

tree16=./INPUT_DATA/NR99_otus_16S.tre
tree18=./INPUT_DATA/tree18.tre

gradient='fresh,brackish,marine'
gradientNAME=SalinityEnviron

GradBinFP=./qtag_code/QTAG.py
RscriptFP=./qtag_code/QTAG_graphing.R

getOTUs=./downstream_rcode/02_make_list_sharedBF_OTUs.R
basicmb=./downstream_rcode/03a_Basic_MB_stats.R
tolerRange=./downstream_rcode/03b_tolerance_ranges.R
taxasum=./downstream_rcode/03c_TaxonomicSummaries_all.R
typetaxa=./downstream_rcode/03d_TypesAcrossTaxonomy.R
salpd=./downstream_rcode/03e_PhyloDistance.R
congruencylvls=./downstream_rcode/03f_Congruency_through_levels.R
brackotus=./downstream_rcode/03g_BrackishOTUs.R

################

dobaltic='True'
dofraser16='True'
dofraser18='True'
doall='True'


mkdir OUTPUT 

if [ $dobaltic == 'True' ] 
then
	# This is for Baltic 16s

	# Make sure is \n
	tr '\r\n' '\n' < $MFBALTIC > temp.txt
	tr '\n\n' '\n' < temp.txt > $MFBALTIC
	rm temp.txt

	python3 $GradBinFP -t $biomBALTIC -m $MFBALTIC -M $gradientNAME --gradient $gradient -o ${outputBALTIC} -R ../$RscriptFP

fi


if [ $dofraser16 == 'True' ]
then

	# This is for Fraser 16s

	# Make sure is \n
	tr '\r\n' '\n' < $MFFRASER16 > temp.txt
	tr '\n\n' '\n' < temp.txt > $MFFRASER16
	rm temp.txt

	python3 $GradBinFP -t $biomFRASER16 -m $MFFRASER16 -M $gradientNAME --gradient $gradient -o ${outputFRASER16} -R ../$RscriptFP

fi


if [ $dofraser18 == 'True' ]
then

# This is for Fraser 18s
	# Make sure is \n
	tr '\r\n' '\n' < $MFFRASER18 > temp.txt
	tr '\n\n' '\n' < temp.txt > $MFFRASER18
	rm temp.txt

	python3 $GradBinFP -t $biomFRASER18 -m $MFFRASER18 -M $gradientNAME --gradient $gradient -o ${outputFRASER18} -R ../$RscriptFP

fi 

## Make a log file with names of folders so I don't have to change it all the time
echo OUTPUT/$outputBALTIC > OUTPUT/dirNames.txt
echo OUTPUT/$outputFRASER16 >> OUTPUT/dirNames.txt
echo OUTPUT/$outputFRASER18 >> OUTPUT/dirNames.txt

mv $outputBALTIC ./OUTPUT/$outputBALTIC
mv $outputFRASER16 ./OUTPUT/$outputFRASER16
mv $outputFRASER18 ./OUTPUT/$outputFRASER18



if [ $doall == 'True' ] 
then
	### Now, get list of otus
	Rscript $getOTUs
	
	### make trees
	
	# For Baltic
	filter_tree.py -i $tree16 -t OUTPUT/${outputBALTIC}/OTUs.txt -o OUTPUT/tree_B16_filt.tre

	# For Fraser
	filter_tree.py -i $tree16 -t OUTPUT/${outputFRASER16}/OTUs.txt -o OUTPUT/tree_F16_filt.tre

	## New tree
	filter_tree.py -i $tree18 -t OUTPUT/${outputFRASER18}/OTUs.txt -o OUTPUT/tree_F18_filt.tre


	# For combo baltic and fraser 16s; **** must already have made non-duplicate OTU list from both datasets
	filter_tree.py -i $tree16 -o OUTPUT/tree_BF16_filt.tre -t OUTPUT/allOTUs_combo_BF16.txt


	# Run remaining scripts
	Rscript $basicmb
	Rscript $tolerRange
	Rscript $taxasum
	Rscript $typetaxa
	Rscript $salpd
	Rscript $congruencylvls
	Rscript $brackotus


fi



