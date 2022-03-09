#barcode_statistics.FDR.Aggregate               execfile("barcode_statistics.FDR.py")
#python code can aggregate outputs from the FDR / FET scripts int a more usable form
#
# 2017.06.15 - Updated to go therouth FET output to get pvalue and table
# 2017.06.06 - Output pvalues -  but - outputs all adj pvalue there is no pvalue in those sheets
# 2017.04.12 - Create -  aggregate all FDR values (padj) into one df, then save

import csv
import pandas as pd
import time
date = time.strftime("%Y.%m.%d")+"_"

#--------------------------------------
# open file with list of all files
with open('FET_files.txt', 'rt') as fin: # this is just a list of the file names that will go through
        cin = csv.reader(fin, delimiter = ',')
        FET_files = [row for row in cin]

# df to store all data
#df_all_adj_pvalues = pd.DataFrame()
#df_all_adj_pvalues_sum_by_gene = pd.DataFrame()
#df_all_adj_pvalues_count_by_condition = pd.DataFrame()
df_all_pvalues = pd.DataFrame()
df_all_table = pd.DataFrame()

#Loop through all files +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for FET_file in FET_files:

    # open file to do FET on
    df_FET_file = pd.read_csv(FET_file[0], skipinitialspace=True, header=range(7), index_col=0)

    # concat padj data with main file
    #df_all_adj_pvalues = pd.concat([df_all_adj_pvalues, df_FET_file.filter(regex='adj')], axis=1)

    # concat p data with main file #
    #df_all_pvalues = pd.concat([df_all_pvalues, df_FET_file.filter(regex='pvalue')], axis=1)
    #df_all_pvalues = pd.concat([df_all_pvalues, df_FET_file.xs(['4','50','pvalue'],level=[1,4,5],axis=1,drop_level=False)], axis=1)

    # concat padj sum data with main file
    #df_all_adj_pvalues_sum_by_gene = pd.concat([df_all_adj_pvalues_sum_by_gene, df_FET_file.filter(regex='sum')], axis=1)

    # concat table data with main file
    df_all_table = pd.concat([df_all_table, df_FET_file.xs(['4','50','table'],level=[1,4,5],axis=1,drop_level=False)], axis=1)

    # aggregate data for each condition
    #df_all_adj_pvalues_count_by_condition = df_all_adj_pvalues_count_by_condition.append(pd.melt(df_FET_file.filter(regex='adj_').ix[0:1]), ignore_index=True)


    print FET_file[0], df_FET_file.filter(regex='adj_').shape

# save files
#df_all_adj_pvalues.to_csv(FET_file[0][:11]+'all_adj_pvalues.csv')
#df_all_pvalues.to_csv(FET_file[0][:11]+'all_pvalues.csv')
#df_all_adj_pvalues_sum_by_gene.to_csv(FET_file[0][:11]+'all_adj_pvalues_sum_by_gene.csv')
df_all_table.to_csv(FET_file[0][:11]+'all_tables.csv')
#df_all_adj_pvalues_count_by_condition.to_csv(FET_file[0][:11]+'all_adj_pvalues_count_by_condition.csv')
