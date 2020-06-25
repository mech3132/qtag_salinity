#!bin/bash

### To re-create json versions
### ONLY RUN IF ON REMOTE SERVERS FORWHICh HDF5 DOESN'T WORK
### 16S

grep "Chloroplast" Fraser_16S/OTUTable_text.txt > Fraser_16S/toDelete.txt
grep "Mitochond" Fraser_16S/OTUTable_text.txt >> Fraser_16S/toDelete.txt
grep "Eukary" Fraser_16S/OTUTable_text.txt >> Fraser_16S/toDelete.txt

## Manual removal of contaminants; contaminants from contam analysis. ##
grep "GQ350231.1.1372" Fraser_16S/OTUTable_text.txt >> Fraser_16S/toDelete.txt
grep "JX525782.1.1434" Fraser_16S/OTUTable_text.txt >> Fraser_16S/toDelete.txt
grep "KY608105.1.1200" Fraser_16S/OTUTable_text.txt >> Fraser_16S/toDelete.txt

biom convert -i Fraser_16S/OTUTable_text.txt --to-json --table-type="OTU table" --header-key taxonomy -o Fraser_16S/OTUTable_16S.biom 

rm Fraser_16S/otu_table_nochloromito_wtaxa.biom
filter_otus_from_otu_table.py -i Fraser_16S/OTUTable_16S.biom -e Fraser_16S/toDelete.txt -o Fraser_16S/otu_table_nochloromito_wtaxa.biom

# Now, let's collapse samples
rm Fraser_16S/MF_16sFraser_noConCOL.txt
rm Fraser_16S/otu_table_nochloromito_col_wtaxa.biom
collapse_samples.py -b Fraser_16S/otu_table_nochloromito_wtaxa.biom -m Fraser_16S/Fraser16s_mappingfile_merged.txt --output_biom_fp Fraser_16S/otu_table_nochloromito_col_wtaxa.biom --output_mapping_fp Fraser_16S/MF_16sFraser_noConCOL.txt --collapse_fields 'ColRep'


### 18S

rm Fraser_18S/otu_table_col_wtaxa.biom
biom convert -i Fraser_18S/OTUTable_text.txt --table-type="OTU table" --to-json --header-key taxonomy -o Fraser_18S/otu_table_col_wtaxa.biom

