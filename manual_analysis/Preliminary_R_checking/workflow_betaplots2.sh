#!/bin/bash

##### Preliminary beta div plots to check where all samples lay ####
# Ensure that you are in the macqiime environment already

## Paths for MC Computer
fraser16=../../INPUT_DATA/Fraser_16S/otu_table_nochloromito_wtaxa.biom
fraser18=../../INPUT_DATA/Fraser_18S/otu_table_nocon.biom
tree16=../../INPUT_DATA/NR99_otus_16S.tre
tree18=../../INPUT_DATA/RAxML_labelledTree.EPA_placement_18s.figtree.tre
followupScript=collapsedstats.sh
plotscript=PCOA_plots.R


### START ####

sed -e 's/QUERY___//g' $tree18 > ./tree18.tre

filter_samples_from_otu_table.py -i $fraser16 -n 1000 -o fraser16_otutable_min1000.biom
filter_samples_from_otu_table.py -i $fraser18 -n 1000 -o fraser18_otutable_min1000.biom

single_rarefaction.py -i fraser16_otutable_min1000.biom -d 1000 -o fraser16_otutable_rare1000.biom
single_rarefaction.py -i fraser18_otutable_min1000.biom -d 1000 -o fraser18_otutable_rare1000.biom


beta_diversity.py -i fraser16_otutable_min1000.biom -t $tree16 -m bray_curtis,weighted_unifrac,unweighted_unifrac -o beta_div_fr16_min
beta_diversity.py -i fraser16_otutable_rare1000.biom -t $tree16 -m bray_curtis,weighted_unifrac,unweighted_unifrac -o beta_div_fr16_rare

beta_diversity.py -i fraser18_otutable_min1000.biom -t tree18.tre -m bray_curtis,weighted_unifrac,unweighted_unifrac -o beta_div_fr18_min
beta_diversity.py -i fraser18_otutable_rare1000.biom -t tree18.tre -m bray_curtis,weighted_unifrac,unweighted_unifrac -o beta_div_fr18_rare

sh $followupScript

Rscript $plotscript

