#!bin/bash sh

## Paths for MC Computer
otu16=../../INPUT_DATA/Fraser_16S/otu_table_nochloromito_col_wtaxa.biom
otu18=../../INPUT_DATA/Fraser_18S/otu_table_col_wtaxa.biom
mf16=../../INPUT_DATA/Fraser_16S/MF_16sFraser_noConCOL.txt
mf18=../../INPUT_DATA/Fraser_18S/MF_18sFraser_noConCOL.txt
tree16=../../INPUT_DATA/NR99_otus_16S.tre
tree18=tree18.tre

single_rarefaction.py -i $otu16 -o fraser16_OTUTable_rare1000_col.biom -d 1000
single_rarefaction.py -i $otu18 -o fraser18_OTUTable_rare1000_col.biom -d 1000


beta_diversity.py -i fraser16_OTUTable_rare1000_col.biom -m weighted_unifrac,unweighted_unifrac,bray_curtis -t $tree16 -o beta_div_COLLAPSED_16S
beta_diversity.py -i fraser18_OTUTable_rare1000_col.biom -m weighted_unifrac,unweighted_unifrac,bray_curtis -t $tree18 -o beta_div_COLLAPSED_18S