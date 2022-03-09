#barcode_statistics.FDR		execfile("barcode_statistics.FDR.py")
#python code to do false discovery rate correction on a set of pvalues
#
# 2017.04.12 - Complete re-write to use Pandas and to drop genes with different number of alleles
# 2016.03.18b - Output into one file with adjusted p-value in column to right of p-value
# 2016.03.18 - File input, working
# 2016.03.17 - Create

import csv
import numpy as np
import pandas as pd
import ast

# for R package to work
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
R_stats = importr('stats')


def FDR_adjust_pvalues(pvalue_list, N=None, method='BH'):
    """ Adjust a list of p-values for false discovery rate using R's stats::p.adjust function.

    N and method are passed to R_stats.p_adjust:
     - N is the number of comparisons (if left unspecified, defaults to len(pvalue_list), I think)
     - method is the name of the adjustment method to use (inherited from R)

    Note that this MUST be done after all the p-values are already collected, on the full list of p-values at once:
     trying to do it on single p-values, even with adjusted N, will give different results!
    """
    if not method in R_stats.p_adjust_methods:
        raise ValueError("Unknown method %s - method must be one of (%s)!"%(method, ', '.join(R_stats.p_adjust_methods)))
    if N is None:   return R_stats.p_adjust(FloatVector(pvalue_list), method=method)
    else:           return R_stats.p_adjust(FloatVector(pvalue_list), method=method, n=N)

def table_count(x,num_allels):
    if sum(ast.literal_eval(x.replace(":",","))) > num_allels:
        return True
    else:
        return False

def tup_replace(tup, val, ix): # replace tublpe value at specific index location
    return tup[:ix] + (val,) + tup[ix+1:]

#--------------------------------------
# Variables
adj_pvalue_threshold = 0.3

#--------------------------------------
# open file with list of all files
with open('FET_files.txt', 'rt') as fin: # this is just a list of the file names that will go through
	cin = csv.reader(fin, delimiter = ',')
	FET_files = [row for row in cin]

#Loop through all files +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for FET_file in FET_files:
    
    # DO FDR ON ALL P-VALUES, using allel cuttoff
    
    # open file to do FET on
    df_FET_file = pd.read_csv(FET_file[0], skipinitialspace=True, header=range(7), index_col=0)
    
    #cols to iterate through
    cols_table = df_FET_file.filter(regex='table').columns
    cols_pvalue = df_FET_file.filter(regex='pvalue').columns

    #make output df
    df_FET_file_and_FDR = df_FET_file

    # loop through number of allels needed to preform FDR (>) #-----------------------------------------------
    for num_allels in [0,1,2,3]:

        # create bool df if gene meets IB threshold
        df_bool = df_FET_file.filter(regex='table').replace(to_replace=":",value=",",regex=True).fillna('[0]').applymap(lambda x: table_count(x,num_allels)) # old way .applymap(table_count(x,num_allels))

        # creat df to store adj pvalues
        df_FDR_pvalues = pd.DataFrame(index =  df_FET_file.index, columns = cols_pvalue)
        columns = list(cols_pvalue.values) # need to do this to disconect columns from the cols_pvalue, otherwise when changing it later it will change!
        df_FDR_pvalues.columns=pd.MultiIndex.from_tuples(columns)
        
        # create list to store # genes in each correction
        l_FDR_gene_count = []
        
        # loop through all pvalues doing FDR correction
        for col_pvalue, col_table in zip(cols_pvalue, cols_table):
            # make sure col in pvalue and table the same!
            if col_pvalue[:5] <> col_table[:5]:
                print 'not the same columns!'
                break
            # do fdr on pvalue col that meet number of allels
            l_FDR_pvalues = list(FDR_adjust_pvalues(df_FET_file.ix[df_bool[col_table],col_pvalue])) # dropping na values does not effect FDR value

            # map the pvalues back to the df (note not all will be there because they were dropped if didn't meet number of allels)
            #     # this makes a dictionary of the pvalue and the index (could just make a series), then updates the column
            #     # based on the values in the Series (dict), values not there are just np.nan
            df_FDR_pvalues[col_pvalue].update(pd.Series(dict(zip(df_FET_file.ix[df_bool[col_table]].index,l_FDR_pvalues))))
            
            # count number of genes less than
            l_FDR_gene_count.append((df_FDR_pvalues[col_pvalue] < adj_pvalue_threshold).sum())
            
        # rename col for adj_pvalue and number of allels
        df_FDR_pvalues.columns.set_levels(['adj_pvalue-alleles_'+str(num_allels)]*len(cols_pvalue), 5, inplace=True)
        
        # add l_FDR_gene_count to header in line 6
        df_FDR_pvalues.columns.set_levels(l_FDR_gene_count, 6, inplace=True)
        df_FDR_pvalues.columns.set_labels(range(len(cols_pvalue)), 6, inplace=True) # change the label for the individual counts to match the level options
        
        # add adj pvalue to orginal df
        df_FET_file_and_FDR = pd.concat([df_FET_file_and_FDR, df_FDR_pvalues], axis=1)

        # sort
        df_FET_file_and_FDR.sortlevel(0, axis=1, inplace = True)
    #------------------------ all alleles ---------------------------------

    # add in last column sum of all hits
    sum_col_name = tup_replace(tup_replace(df_FET_file_and_FDR.columns[0],'sum',5),'sum',6)
    df_FET_file_and_FDR[sum_col_name] = (df_FET_file_and_FDR.filter(regex='adj')< adj_pvalue_threshold).sum(axis=1)

    df_FET_file_and_FDR.to_csv(FET_file[0][:-3] + 'FDR.csv', ) # columns=False, index=False
    
    print 'done with ' + FET_file[0]

print 'done'
