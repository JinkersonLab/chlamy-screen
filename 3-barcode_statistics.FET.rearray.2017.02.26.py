#barcode_statistics		execfile("barcode_statistics.rearray.2017.02.26.py")
#python code to do statistical analysis on barcode frequency 

#The code prepares the ratios then imports them into a function that looks at statistics
#
# **** things to do to get running on server, delete # from import scipy and fischer test lines, make sure sample name ' Plates % etc' is correct
#
#2017.02.26	Make it so that it runs with merged data; add float function to handle '' as 0, and float for confidence
#2016.09.19	Allow input arguemtns to be specified that are the IB file and also samples to compare
#2016.09.06	Specify that Fisher exact test is greater (not defalut two tailed) so will not give probability to under represented genes
#2016.09.06	Fix so that it is checking raw read counts, not normalized read counts
#2016.09.06	do not calculate fisher excat test everytime, check dictionany if it is there, then procede to calculuate if not already calculated
#2016.03.17	included filtering based on: confidence, feature, threshold, and reads. all filtering can be input by samples_to_compare.tab. ONLY ONE THRESHOLD PER LINE (ie make two lines one for > one for <)
#2016.03.14	- work on getting output right. need to put in threshold, read count, and feature type limits
#2015.11.25 	working on server w Fisher test, need to make output files, etc, multiple threholds and check to make sure it makes sense
#2015.11.24 	intial writing

#INPUTS
#	1) spreadsheet with data
#	2) samples_to_compare.tab -> tab spreadsheet with sample in column 1 and control (divided by) in column 2

#OUTPUTS
#
#---headers--
#sample_v_sample		:samples that were compaired (first divided by second)
#confidence <=		:the mutant confidence number threshold that a barcode had to meet to be included in the analysis. 5 = all mutants; 3 = less than or equal to 3 needed
#feature				:feature barcode/mutant had to be to be included (all = any feature) ['CDS', '5UTR'] =  only these features
#threshold			:threshold that sample/control must meet to be a 'hit'
#min_num_reads >=		:minimum number of reads needed to include barcode in analysis
#column				:table, odds_ratio, or pvalue
#all_genes			:all the genes that were included in the analysis for given criteria [hits, not hits]

#----------------------------------------------
# Initialize
#----------------------------------------------
import csv
import time
import subprocess
import operator
import ast	#for converting strings that look like list into real lists
import numpy as np
from collections import defaultdict
import scipy.stats#E!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!undo to run on server!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
import sys

def f_check_feature(f_feature, f_feature_type): #input: 1) feature from row, 2)features to check #this function checks to see if the feature that is in the given row is of the feature type criteria
	if f_feature_type == 'all':
		return True
	else:
		for feature_type in f_feature_type:
			if feature_type == f_feature.strip():
				return True
		#didn't find one, so return false
		return False	

def f_threshold_operator(f_sample, f_control, f_operator, f_threshold): #function to determine if threshold has been met. Returns True or False
	if f_operator == '>':
		return f_sample / f_control > f_threshold
	elif f_operator == '>=':
		return f_sample / f_control >= f_threshold
	elif f_operator == '<':
		return f_sample / f_control < f_threshold
	elif f_operator == '<=':
		return f_sample / f_control <= f_threshold

def mk_float(s): # use this isntad of float, beacause if there is '' then it will give you 0, where float('') gives error!
    s = s.strip()
    return float(s) if s else 0
	
# Constants-----------------------------------------------
gene_name_column = 4 # The column that contains the gene name
gene_feature_column = 10 # The column that contains the gene feature 'gene name; feature; confidence'
gene_confidence_column = 11 # this is from rearrayed data where the column is different

# Default inputs -----------------------------------------> commented out so that in file inputs used
#thresholds = [5,10] # fold threshold
#reads = [10,500] # calculate statistics if reads are more than minimum
#l_feature = ['all',['CDS','5UTR']] # which features to allow: intergenic, intron, 3UTR, CDS, 5UTR, MULTIPLE_SPLICE_VARIANTS,
#l_confidence = [5, 4, 3] # confidence values less than or equal to these pass through
#gt_or_lt_operator = '>' # change this for changing if you are looking at less than or greater than thresholds
gt_or_lt_operator_name = ['na'] # this initializes the operator name

# Variables------------------------------------------------
all_IB = []
d_genes ={}

date = time.strftime("%Y.%m.%d")+"_"

#----------------------------------------------------------
# Load files
#----------------------------------------------------------
# Look for user supplied arguemtns, if not then use defaults
if len(sys.argv) > 1: # does file have arguments?
	try:
		if sys.argv[1][-3:] == 'csv':
			IB_data_file = sys.argv[1]
			samples_to_compare_file = sys.argv[2]
		elif sys.argv[1][-3:] == 'tab':
			samples_to_compare_file = sys.argv[1]
			IB_data_file = sys.argv[2]
	except: 
		print "not enough input arguments, 2 required: IB data file and samples to compare file"
else:
	IB_data_file = 'R123456.merged-5p3p.max-initial.best-IB-conf-feat-2017.02.26.csv'
	samples_to_compare_file = 'samples_to_compare.lt.tab'

#Load output data (normalized read counts)
with open(IB_data_file, 'rt') as fin:# on computer: 2015.11.17_all.IB.test.csv   on server: 2016.03.15_all.IB.csv    \\\\Full len on computer 2015.11.17_all.IB.csv ///'2015.11.17_all.IB.csv' #server file: 2015.11.18_all.IB.csv !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
	cin = csv.reader(fin, delimiter = ',')
	all_IB = [row for row in cin]

#load all gene names into dict (for output)
for row in all_IB:
	d_genes[row[gene_name_column]] = []

#try to import to numpy right away
#all_IB2 = np.loadtxt('2015.11.17_all.IB.csv', delimiter=',')

#Load file with comparisons to be made
#		*File has 1st two columns with each row having the comparison, column 1 divided by column 2 (sample by control)
#		*File has following up columns that have the following columns
#			*Sample	Control	op	Thresholds	Reads	Features	mut. confidence: example: R2-Anoxia	R2-Initial_IV	<	[0.3,0.1,0.05,0.025,0.01]	[1,5,10,50,100]	[all,['CDS','5UTR']]	[5,4,3]
#		*Inputs should like like lists [1, 2, 3]
with open(samples_to_compare_file, 'rt') as fin:
	cin = csv.reader(fin, delimiter = '\t')
	samples_to_compare = [row for row in cin]

#v################################################### LOOP ########################################3
for samples in samples_to_compare:
	#varable to initalize
	d_header ={'sample_v_sample':['sample_v_sample'], 'threshold':['threshold'], 'column':['column'], 'min_num_reads':['min_num_reads >='], 'feature':['feature'], 'confidence':['confidence <='], 'all_genes':['all_genes']}
	
	#load parameters in file: load op, threshold, reads, features, mut. confidence here (for each sample) sample = [sample, control, op, thresholds, reads, features, mutant confidence] *should look like lists with *[x, y, z]
	gt_or_lt_operator = samples[2] # op
	thresholds = ast.literal_eval(samples[3]) # thresholds (turn string that looks like list into list)
	reads = ast.literal_eval(samples[4]) # reads
	l_feature = ast.literal_eval(samples[5])# feature
	l_confidence = ast.literal_eval(samples[6])# confidence
	
	# Find the columns that have the samples to compare
	sample_expt_loc = all_IB[0].index(samples[0] +'_normalized_reads') # Column 1 = expt (have to add on +' Plates % of all reads' so can find normalized data)
	sample_control_loc = all_IB[0].index(samples[1] +'_normalized_reads') # Column 2 = control ### computer version ' % of all reads'
	sample_control_reads_loc = all_IB[0].index(samples[1] +'_read_count') # Column where the raw reads are found
	
	for confidence in l_confidence: # go through confidence thresholds
		for feature in l_feature: #	go through the features
			for min_threshold in thresholds:	
				for min_num_reads in reads:	
					#Varibales to be new for every sample
					d_gene_ratios = defaultdict(list)
					d_gene_fisher = {} # dictionary with [0,0] list where 0 = # of genes that meets cutoff; 1 = # of genes that does not meet cutoff
					d_all_fisher = [0,0]
					
					#----------------------------------------------
					# Prepare ratios
					#----------------------------------------------
					
					#Find the columns that have the samples to compare was here....
					#sample_expt_loc = all_IB[0].index(samples[0] +' Plates % of all reads') # Column 1 = expt (have to add on +' Plates % of all reads' so can find normalized data)
					#sample_control_loc = all_IB[0].index(samples[1] +' Plates % of all reads') # Column 2 = control
					
					#Divide columns by each other and load into dictionary (divide 0 by 1 (column))
					for row in all_IB[1:]:
						if float(row[gene_confidence_column]) <= confidence: # check that IB meets confidence threshold (confidence is 3rd item in column)
							if f_check_feature(row[gene_feature_column], feature): # check to make sure feature is approved- function input 1) feature from row (split and reomve white space)	 2) features allowed ('all' or ['CDS', '5UTR'], etc)
								if mk_float(row[sample_control_reads_loc]) >= min_num_reads:# if more reads than given values, then(was <>0)If control is 0 then do not process here		float(row[sample_0_loc]) <> 0 and float(row[sample_1_loc]) <> 0
									d_gene_ratios[row[gene_name_column]].append(mk_float(row[sample_expt_loc])/mk_float(row[sample_control_loc]))
									
									#fisher exact test prep
									d_gene_fisher.setdefault(row[gene_name_column],[0,0]) #set default for gene as [0, 0]
									
									if f_threshold_operator(mk_float(row[sample_expt_loc]), mk_float(row[sample_control_loc]), gt_or_lt_operator, min_threshold): #float(row[sample_expt_loc])/float(row[sample_control_loc]) >= min_threshold: #add genes to column 0 if greater than thershold
										d_gene_fisher[row[gene_name_column]][0] += 1
										d_all_fisher[0] += 1 # all of the genes
									else:
										d_gene_fisher[row[gene_name_column]][1] += 1 # add genes to second column if not meeting threshold
										d_all_fisher[1] += 1 # all of the genes
										
								elif mk_float(row[sample_expt_loc]) <> 0 and mk_float(row[sample_control_loc]) == 0:
									pass
									#Infinity (no reads in control, only expt)
									#d_gene_ratios[row[gene_name_column]].append('INF')
								
					#----------------------------------------------
					# Function to calculate statistics
					#----------------------------------------------
					#
					#		Input	- dictonary with gene name as key and list of all ratios
					#				- all values
					
					# Create Variables
					d_gene_fisher_pvalue = {}
					d_threshold_list = {}
					
					# Fisher exact test
					
					for gene, threshold_list in d_gene_fisher.items():
						d_gene_fisher_pvalue.setdefault(gene,[0,0])
						
						if tuple(threshold_list) not in d_threshold_list:
							d_threshold_list[tuple(threshold_list)] = scipy.stats.fisher_exact([threshold_list, d_all_fisher], alternative='greater') #alternative='two-sided'  gives you p values for both sides
						
						# Combine gene info wit Fisher Exact test values
						d_gene_fisher_pvalue[gene][0], d_gene_fisher_pvalue[gene][1] = d_threshold_list[tuple(threshold_list)][0],d_threshold_list[tuple(threshold_list)][1]
						
					#---------------------------------------
					# OUTPUT
					#---------------------------------------
					
					#build headers (need 3 X of everything) (every once thorugh generates threholds, odds ratio, and pvalue)
					d_header['column'].extend(['table','odds_ratio','pvalue'])
					d_header['threshold'].extend([gt_or_lt_operator + str(min_threshold), gt_or_lt_operator + str(min_threshold), gt_or_lt_operator + str(min_threshold)])
					d_header['min_num_reads'].extend([min_num_reads, min_num_reads, min_num_reads])
					d_header['feature'].extend([feature, feature, feature])
					d_header['confidence'].extend([confidence, confidence, confidence])
					d_header['sample_v_sample'].extend([(samples[0]+' v '+samples[1])]*3)
					
					#output all genes (in header dictionary)
					d_header['all_genes'].extend([d_all_fisher, 'NA', 'NA'])
					
					#put stats data into dictionary with all gene names (need to include all gene names and then put na if no data available, so that if next run may be data will be aligned)
					#3 values per time [thresholds], odds ratio, pvalue (need to have 3 times of everything for header)
					for gene, data in d_genes.items():
						try:
							data.extend([d_gene_fisher[gene], d_gene_fisher_pvalue[gene][0], d_gene_fisher_pvalue[gene][1]])
						except:
							data.extend(['NA','NA','NA'])
			
	#determine > or < for outfile
	if gt_or_lt_operator == '>':
		gt_or_lt_operator_name = 'gt'
	elif gt_or_lt_operator == '>=':
		gt_or_lt_operator_name = 'gteq'
	elif gt_or_lt_operator == '<':
		gt_or_lt_operator_name = 'lt'
	elif gt_or_lt_operator == '<=':
		gt_or_lt_operator_name = 'lteq'

	
	# Output headers and data
	with open(date+(samples[0]+'_v_'+samples[1])+ '.' + gt_or_lt_operator_name + '.FET.csv', 'wb') as fout: #was 'wt' *for windows 'wb' will not include extra returns in csv file (wt does)
		csv_writer = csv.writer(fout)
		for header in ['sample_v_sample', 'confidence', 'feature', 'threshold', 'min_num_reads', 'column', 'all_genes']: #adds these headers in this order at top of file
			csv_writer.writerow(d_header[header]) #adds header to top of file
		for key, value in d_genes.items(): # goes through dictionary and writes each row
			value.insert(0,key)
			csv_writer.writerow(value)
	# Reset for next sample comparison 
	
	d_genes = dict((k,list()) for k in d_genes) #resets the dictionary with empty lists
	print 'done with '+samples[0]+'_v_'+samples[1]
#execfile("barcode_statistics.py")

print 'done'
