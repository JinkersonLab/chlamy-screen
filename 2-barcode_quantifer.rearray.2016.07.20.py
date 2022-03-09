#Barcode Quantifier   execfile("barcode_quantifer_rearray_2016.07.20.py")
#
#Inputs
#1)Takes barcode quantification output
#2)Spread sheet detail file with barcode and plate/well locations
#
# IF HAVING ARRAY PROBLEM! Having an numpy array problem, check to make sure barcode (5prime or 3prime) is correct!!!!
#
#2016.07.20 Update it to split the gene name, confidence, feature type
#2016.07.18 Update for the rearray library file (large-lib_rearray_2016.07.10.txt), have to update the plates becasue rearray plates have C in name
#2015.09.30	Can loop through files and draw figs.
#2015.09.15	Can loop (need input file: index_name_description.txt); fixed divide by zero error
#2015.09.03	Works with test input, output, added % and headers to output
#
#Bugs
#output has 1 then 2 then 3 samples (they are added sequentally, need to take last set which has all of the samples)
#
#
import csv
import string
import numpy as np
import matplotlib
matplotlib.use('pdf') #TkAgg # Use 'Agg' for png (can use multiple .use() invocations to use more than one)
matplotlib.use('agg')
import matplotlib.pyplot as plt
import subprocess
from collections import defaultdict #to make dict with list is defalut
import time

#Variablesf
start_plate_numb = 8
end_plate_numb = 245
plates = [["C"+str(x).zfill(3),0,0,0] for x in range(start_plate_numb,end_plate_numb+1)] #oldway, ["plate"+str(x),0,0,0]
all_plates = [["C"+str(x).zfill(3)] for x in range(start_plate_numb,end_plate_numb+1)] # oldway ["plate"+str(x)]
barcode_end = "3prime" # or "3prime" 5prime
date = time.strftime("%Y.%m.%d")+"_"

number_of_barcodes = 0
sum_of_reads = 0

wells = []
all_wells = []
#make list of wells:wq

for letter in string.ascii_uppercase[0:16]: # make wells
	for num in ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24']: #range(1,25):
		wells.append([letter+str(num),0,0])
		all_wells.append([letter+str(num)])

all_plates_header = [['plates'],['name']]
all_wells_header = [['wells'],['name']]

class Counter(dict):
	def __missing__(self, key):
		return 0
#-----------------------------------------------------------------------------
# Load files (once only)
# 	import file with names and description, tab seperated two coloums, index [0], description per row[1], description two[2], plate start(for analysis)[3], plate end (for analysis)[4]
with open('index_name_description.txt', 'rt') as fin:
	cin = csv.reader(fin, delimiter = '\t')
	index_name_description = [row for row in cin]
#determine number of samples
numb_samples = len(index_name_description) #-1 if only one sample

#Load library data from Weronika (Barcode 17, plate 0, well 1, barcode end 2, cre 7, gene name 22) #old file large-lib_fixed-positions_RJ.txt
with open('large-lib_rearray_2016.07.10_RJ.txt', 'rt') as fin:
	cin = csv.reader(fin, delimiter = '\t')
	library_data = [row for row in cin]
#-----------------------------------------------------------------------------
	
#Makes dictionary from  Weronika Data
d_barcode = {}
d_barcode_count_sum = {} #put all barcodes into one dictionary
d_plate_well_gene = defaultdict(list) #Put all the gene location for 5p and 3p at well location; will assign at end
d_plate_well_gene_condensed = defaultdict(list) # puts all the gene location for 5 and 3p

#Make Counters only once
libary_barcodes_per_plate = Counter()
libary_barcodes_per_well = Counter()

#Count barcodes in library on each plate/well
for row in library_data:
	#d_plate_well_gene[row[1]+row[2]].append([row[10],row[35]+"; "+row[12]+"; "+row[16],row[5]]) #cre, gene name; feature, 5/3prime #WERONIKA LIBRARY COLUMN; #d_plate_well_gene[row[0]+row[1]].append([row[7],row[26]+"; "+row[9]+"; "+row[13],row[2]]) #ORIGINAL
	d_plate_well_gene[row[1]+row[2]].append([row[10],row[34], row[35], row[36],row[44],row[46], row[12], row[16],row[5]]) #cre, gene name; feature, 5/3prime #WERONIKA LIBRARY COLUMN; #d_plate_well_gene[row[0]+row[1]].append([row[7],row[26]+"; "+row[9]+"; "+row[13],row[2]]) #ORIGINAL
	if row[5] == barcode_end: # only add if it matches barcode end (5prime or 3prime) 
		# IB : plate; well; side; gene; gene name;
		d_barcode.update({row[29]:[row[1],row[2],row[5],row[10],row[34], row[35], row[36],row[44],row[46], row[12], row[16],[[0,0]]*numb_samples]}) #dictonary with barcode as key and locations etc as value; last list is where data goes #original d_barcode.update({row[21]:[row[0],row[1],row[2],row[7],row[26]+"; "+row[9]+"; "+row[13],[[0,0]]*numb_samples]}) #dictonary with barcode as key and locations etc as value; last list is where data goes
		#sum up barcodes on each plate/well in library file
		libary_barcodes_per_plate[row[1]] += 1
		libary_barcodes_per_well[row[2]] += 1


l_d_barcode_header = [library_data[0][29],library_data[0][1],library_data[0][2],library_data[0][5],library_data[0][10],library_data[0][34], library_data[0][35], library_data[0][36], library_data[0][44], library_data[0][46], library_data[0][12], library_data[0][16]]

#Make Summary list
l_summary = [["sample"],["# of barcodes in sequencing"],["# of barcodes found in library"],["# of barcodes found in library in select plates" ], ["sum of reads from all barcodes"],["sample reads"],["sum of reads from barcodes found in library"],["sum of reads from barcoes found in select plates"]]

#-----------------------------------

#make figures


#----- fig_hist
fig_len = numb_samples * 2.5 #11 inches

fig_hist, ax_hist = plt.subplots(numb_samples,2) #makes figures and array with plots ( # of rows, # of columns) START WITH 1 NOT 0
fig_hist.subplots_adjust(left=0.1, right=0.9,
	bottom=0.1, top=0.9,
	hspace=0.4, wspace=0.4) #hspace = height reserved for white space between subplots; the amount of width reserved for blank space between subplots
fig_hist.set_size_inches(8.5,fig_len)

#----- fig_plate
fig_plate, ax_plate = plt.subplots(numb_samples,3) #(rows, columns)
fig_plate.subplots_adjust(left=0.1, right=0.9,
	bottom=0.1, top=0.9,
	hspace=0.4, wspace=0.4)
fig_plate.set_size_inches(8.5,fig_len)

#----- fig_wells
fig_wells, ax_wells = plt.subplots(numb_samples,2) #(rows, columns)
fig_wells.subplots_adjust(left=0.1, right=0.9,
	bottom=0.1, top=0.9,
	hspace=0.4, wspace=0.4)
fig_wells.set_size_inches(8.5,fig_len)


#-----------------------------------------------------------------------------------------------------------------------------------------------------------
# Loop through all of the samples
for idx, sample in enumerate(index_name_description): #need [:-1] for 1 sample
	
	print sample
	#Counters
	counter_plates = Counter()
	counter_well = Counter()
	counter_cre = Counter()
	counter_end = Counter()
	counter_desc = Counter()
	
	reads_plates = Counter()
	reads_well = Counter()
	reads_cre = Counter()
	reads_end = Counter()
	reads_desc = Counter()
	
	#Lists
	l_barcode_count_found=[]
	l_barcode_count_found_in_plates=[]
	
	#variables
	plate_start_for_analysis = sample[3] #refers to row in index_name_description ### for all plates this should be 101
	plate_end_for_analysis = sample[4] #refers to row in index_name_description ### for all plates this should be 650
	
	#Loads files (multiple times)
	# Loads Barcode count file
	with open(sample[0] + '.count.sort', 'rt') as fin:
			cin = csv.reader(fin, delimiter = ' ')
			l_barcode_count = [[row[0],int(row[1])] for row in cin] #make list
	
	#Sum up total number of reads and barcodes in .count.sort file
	for row in l_barcode_count:
		number_of_barcodes += 1 #count the number of barcodes
		sum_of_reads += int(row[1])
	
	#determine # reads
	np_l_barcode_count = np.array(l_barcode_count, dtype=object)
	sample_reads = sum(np_l_barcode_count[:,1])
	
	# Assign barcode to plate / well; puts count data into d_barcode
	for row in l_barcode_count:
		row.extend(d_barcode.get(row[0], ['NF','NF','NF', 'NF', 'NF', 'NF'])) #checks in dictionary for barcode, adds NF if not found
		if row[2] <> 'NF':
			l_barcode_count_found.append(row) #found in library file(has well location)
			if int(plate_start_for_analysis) <= int(row[2][-3:]) <= int(plate_end_for_analysis):
				l_barcode_count_found_in_plates.append(row) #found in library file and in the plate set for analysis
				d_barcode[row[0]][-1][idx] = [row[1],row[1]*100000000/sample_reads] #[read count; % of reads]# adds to list at end of dic barcode count for sample
	#Make header for d_barcode
	l_d_barcode_header.extend([sample[0]+"_read_count",sample[0]+"_normalized_reads"]) #old way l_d_barcode_header.extend([sample[0]+" "+sample[1]+" read count",sample[0]+" "+sample[1]+" % of all reads"])
	
	# Checks distribution of barcodes
	# Barcodes per plate
	
	for row in l_barcode_count:
		counter_plates[row[2]] += 1
		counter_well[row[3]] +=1
		counter_end[row[4]] += 1
		counter_cre[row[5]] += 1
		counter_desc[row[6]] += 1
		
		reads_plates[row[2]] += int(row[1])
		reads_well[row[3]] += int(row[1])
		reads_end[row[4]] += int(row[1])
		reads_cre[row[5]] += int(row[1])
		reads_desc[row[6]] += int(row[1])
	
	#-----------------
	#Barcodes per plate
	
	#add counts to list of plate
	for row in plates:
		row[-3]=counter_plates[row[0]] # add plate count
		try: #for plates that have 0 barcodes in library, avoid divide by 0 error
			row[-2]=1.0 * counter_plates[row[0]] / libary_barcodes_per_plate[row[0]] # add percent of barcodes found in libary sheet
		except ZeroDivisionError:
			row[-2]=0 # add
		row[-1]=reads_plates[row[0]] # add reads per plate
	
	#makes header for plates output
	plates_header = ["plate", "barcode count", "% of library barcodes found in this plate", "reads per plate"]
	
	#-----------------
	#Barcodes per well
	
	#add counts to list of wells
	for row in wells:
		row[-2] = counter_well[row[0]] # add count of barcodes per well
		#row.extend([1.0 * counter_well[row[0]] / libary_barcodes_per_well[row[0]]]) # add percent of barcodes found in libary sheet
		row[-1] = reads_well[row[0]] # add reads per well
	
	wells_header = ["well", "barcode count", "% of library barcodes found in this well", "reads per well"]
	
	#-----------------
	#OUT Files
	with open(date+sample[0]+'.plates.csv', 'wt') as fout:
		csv_writer = csv.writer(fout)
		csv_writer.writerow(plates_header) # adds header to top of file
		for row in plates: # goes through list and writes each row
			csv_writer.writerow(row)
	
	#OUT Files
	with open(date+sample[0]+'.wells.csv', 'wt') as fout:
		csv_writer = csv.writer(fout)
		csv_writer.writerow(wells_header) # adds header to top of file
		for row in wells:
			csv_writer.writerow(row)
	#--------------------
	
	#Add to master tracker list with all samples
	for row_all_plates, row_plates in zip(all_plates, plates): #master plate list
		row_all_plates.extend(row_plates[1:])
	
	for row_all_wells, row_wells in zip(all_wells, wells): #master well list
		row_all_wells.extend(row_wells[1:])
	
	#Add header to master tracker list
	all_plates_header[0].extend(plates_header[1:])
	all_plates_header[1].extend([sample[0],sample[0],sample[0]])
	
	all_wells_header[0].extend(wells_header[1:]) #add sample[1] = index name description
	all_wells_header[1].extend([sample[0],sample[0],sample[0]])
	
	#PLOT OUTPUT SECTION>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	np_plates = np.array(plates, dtype=object)
	np_wells = np.array(wells)
	
	np_l_barcode_count_found = np.array(l_barcode_count_found, dtype=object)
	np_l_barcode_count_found_in_plates = np.array(l_barcode_count_found_in_plates, dtype=object)
	
	#histogram of #reads / barcode
	hist_lowrange, hist_highrange = 2, 1000
	hist_bins = 50 #hist_highrange / 2
	
	#plt_hist column 0 ----> histogram of # reads / barcodes log; include all barcodes, all found in library, and all found in analysis plates 
	ax_hist[idx, 0].hist(np_l_barcode_count[:,1], bins=hist_bins, log=True)
	ax_hist[idx, 0].hist(np_l_barcode_count_found[:,1], bins=hist_bins, log=True)
	ax_hist[idx, 0].hist(np_l_barcode_count_found_in_plates[:,1], bins=hist_bins, log=True)
	ax_hist[idx, 0].set_xlabel('# of reads')
	ax_hist[idx, 0].set_ylabel('frequency of barcodes')
	ax_hist[idx, 0].set_title(sample[0]+" "+sample[1])
	
	#plt_hist column 1 ----> histogram of # reads / barcodes
	#ax_hist[idx, 1].hist(np_l_barcode_count[:,1], bins=hist_bins,range=(hist_lowrange, hist_highrange)) #has so many hits with 1 read
	ax_hist[idx, 1].hist(np_l_barcode_count_found[:,1], bins=hist_bins,range=(hist_lowrange, hist_highrange))
	ax_hist[idx, 1].hist(np_l_barcode_count_found_in_plates[:,1], bins=hist_bins,range=(hist_lowrange, hist_highrange))
	ax_hist[idx, 1].set_xlabel('# of reads')
	ax_hist[idx, 1].set_title(sample[0]+" "+sample[1])
	#---------------------------------------------------
	
	#fig_plate
	
	#plot #barcodes / plate
	ax_plate[idx,0].plot(np_plates[:,-3]) #plot # barcodes / plate
	ax_plate[idx,0].set_xticks(range(len(np_plates[:,1]))[::100])
	ax_plate[idx,0].set_xticklabels(np_plates[:,0][::100], size='small', rotation='vertical' ) #makes tick marks the name of the well. Skip to included every 24
	ax_plate[idx,0].set_ylabel('# of barcodes')
	
	#plot %barcodes / plate
	ax_plate[idx,1].plot(np_plates[:,-2]) #plot # barcodes / plate
	ax_plate[idx,1].set_xticks(range(len(np_plates[:,1]))[::100])
	ax_plate[idx,1].set_xticklabels(np_plates[:,0][::100], size='small', rotation='vertical' ) #makes tick marks the name of the well. Skip to included every 24
	#ax_plate[idx,0].set_xlabel('Plate #')
	ax_plate[idx,1].set_ylabel('% of barcodes')
	ax_plate[idx,1].set_title(sample[0]+" "+sample[1])
	
	#plot #reads / plate
	ax_plate[idx,2].plot(np_plates[:,-1],'r') #plot # barcodes / plate
	ax_plate[idx,2].set_xticks(range(len(np_plates[:,1]))[::100])
	ax_plate[idx,2].set_xticklabels(np_plates[:,0][::100], size='small', rotation='vertical' ) #makes tick marks the name of the well. Skip to included every 24
	ax_plate[idx,2].set_ylabel('# of reads')
	#--------------------------------------------------
	
	#fig_wells
	
	#plot #barcodes / well
	ax_wells[idx,0].plot(np_wells[:,-2]) #plot # barcodes / well
	ax_wells[idx,0].set_xticks(range(len(np_wells[:,1]))[::24])
	ax_wells[idx,0].set_xticklabels(np_wells[:,0][::24], size='x-small', rotation='vertical') #makes tick marks the name of the well. Skip to included every 
	ax_wells[idx,0].set_ylabel('# of barcodes')
	
	#plot %barcodes / well
	ax_wells[idx,1].plot(np_wells[:,-1]) #plot # barcodes / plate
	ax_wells[idx,1].set_xticks(range(len(np_wells[:,1]))[::24])
	ax_wells[idx,1].set_xticklabels(np_wells[:,0][::24], size='x-small', rotation='vertical') #makes tick marks the name of the well. Skip to included every 
	ax_wells[idx,1].set_xlabel('Well #')
	ax_wells[idx,1].set_ylabel('# reads')
	ax_wells[idx,1].set_title(sample[0]+" "+sample[1])
	
	#NO WELL % of BARCODES CALCULATED
	#plot #reads / plate
	#ax_wells[idx,2].plot(np_wells[:,3],'r') #plot # barcodes / plate
	#ax_wells[idx,2].set_xticks(range(len(np_wells[:,1]))[::100], np_wells[:,0][::100], size='small', rotation=17 ) #makes tick marks the name of the well. Skip to included every 24
	#ax_wells[idx,2].set_ylabel('# of reads')
	
	
	#---- Summary list
	#l_summary = [["sample"],["# of barcodes in sequencing"],["# of barcdoes found in library"],["# of barcodes found in library in select plates" ], ["sum of reads from all barcodes"],["sample reads"],["sum of reads from barcodes found in library"],["sum of reads from barcoes found in select plates"]]
	l_summary[0].append(sample[0]+" "+sample[1])
	l_summary[1].append(number_of_barcodes)
	l_summary[2].append(len(np.array(l_barcode_count_found, dtype=object)[:,1]))
	l_summary[3].append(len(np.array(l_barcode_count_found_in_plates, dtype=object)[:,1]))
	l_summary[4].append(sum_of_reads)
	l_summary[5].append(sample_reads)
	l_summary[6].append(np.array(l_barcode_count_found, dtype=object)[:,1].sum())
	l_summary[7].append(np.array(l_barcode_count_found_in_plates, dtype=object)[:,1].sum())
	
	for row in l_summary:
		print row
	
	#---- RESET for next loop
	number_of_barcodes = 0
	sum_of_reads = 0
	for row in plates:
		row[-1] = 0
		row[-2] = 0
		row[-3] = 0
	for row in wells:
		row[-1] = 0
		row[-2] = 0

# Save figures
fig_plate.tight_layout()
fig_hist.tight_layout()
fig_wells.tight_layout()

fig_hist.savefig(date+'fig_hist.png')
fig_plate.savefig(date+'fig_plate.png')
fig_wells.savefig(date+'fig_wells.png')

fig_hist.savefig(date+'fig_hist.pdf')
fig_plate.savefig(date+'fig_plate.pdf')
fig_wells.savefig(date+'fig_wells.pdf')

#----------------
# All gene locations for plate-well condenceed so every gene name is only there once. Adds the location 5prime or 3prim next to it by concatination
d_temp = {}
for platewell, genes in d_plate_well_gene.items():
	for value in genes:
		if value[0] in d_temp:
			d_temp[value[0]][2] = d_temp[value[0]][2]+ value[2]#in dict already
		else:
			d_temp[value[0]] = value#not in dict already
	genes = []
	value = []
	for key, value in d_temp.items():
		d_plate_well_gene_condensed[platewell].append(value)
	d_temp = {}

#Add condensed gene info to d_barcode list
for key, value in d_barcode.items():
	value.append(d_plate_well_gene_condensed[value[0]+value[1]])

#-----------------
#OUT Files - ALL
with open(date+'all.plates.csv', 'wt') as fout:
	csv_writer = csv.writer(fout)
	csv_writer.writerow(all_plates_header[1]) # adds header to top of file
	csv_writer.writerow(all_plates_header[0]) # adds header to top of file
	for row in all_plates: # goes through list and writes each row
		csv_writer.writerow(row)

with open(date+'all.wells.csv', 'wt') as fout:
	csv_writer = csv.writer(fout)
	csv_writer.writerow(all_wells_header[1]) # adds header to top of file
	csv_writer.writerow(all_wells_header[0]) # adds header to top of file
	for row in all_wells: # goes through list and writes each row
		csv_writer.writerow(row)

# Output of summary file
with open(date+'all.summary.csv', 'wt') as fout:
	csv_writer = csv.writer(fout)
	for row in l_summary:
		csv_writer.writerow(row)

# Out put barcodes
with open(date+'all.IB.csv', 'wt') as fout:
	csv_writer = csv.writer(fout)
	csv_writer.writerow(l_d_barcode_header) #adds header to top of file
	for key, value in d_barcode.items(): # goes through dictionary and writes each row
		csv_writer.writerow([key,value])

# get rid of ' " [ ] 
#"sed -i -e "s/'//g" -e "s/\[//g" -e "s/\]//g" -e "s/\"//g" all.IB.csv"
subprocess.call(["sed","-i","-e","s/'//g","-e","s/\[//g","-e","s/\]//g","-e","s/\"//g",date+"all.IB.csv"])


with open(date+'libary_barcodes_per_plate.csv', 'wt') as fout:
	csv_writer = csv.writer(fout)
	#csv_writer.writerow(l_d_barcode_header) #adds header to top of file
	for key, value in libary_barcodes_per_plate.items(): # goes through dictionary and writes each row
		csv_writer.writerow([key,value])
		
libary_barcodes_per_plate = Counter()
libary_barcodes_per_well = Counter()
