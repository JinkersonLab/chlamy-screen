# chlamy-screen
Code used to quantify barcodes from pooled Chlamydomonas mutant library screens

The following is a set of scripts which can be used to quantify barcodes from a pooled Chlamydomonas mutant library screen. More information on this approach can be found in the following publication: https://www.biorxiv.org/content/10.1101/2020.12.11.420950v1

Scripts:
1-IB_cut_count_sort.sh
This is a bash script that will use cutadapt to find 5' or 3' barcode sequences, count the barcodes, and reverse the barcodes sequences.

2-barcode_quantifer.rearray.py
The barcode quantifier is a python script that takes the barcode quantification output and a spreadsheet (tab delimited) that contains the barcode and plate/well locations as input. It will output summaries of the barcodes detected, including plots.

3-barcode_statistics.FET.rearray.py
This python script prepares the mutant barcode ratios (experiment/control) and then conducts Fisher's exact test on them to determine if any gene has a phenotype. Required inputs are the barcode count data and a samples_to_compare.tab spreadsheet with the samples to compare. Mutants can be included or excluded for the analysis based on confidence, feature type, and minimum number of reads. Thresholds can be set for the ratio (experiment/control) in order to be considered a 'hit'. The Fisher's exact test table and pvalue can be output for each gene for all comparisons made.

4-barcode_statistics.FDR.py
This python script conducts a false discovery rate correction on the set of pvalues generated in the last script. Genes can be excluded from the FDR correction if they did not have enough alleles to be considered.

5-barcode_statistics.Aggregate.py
This python script can aggregate data generated from the previous scripts into more usable forms.
