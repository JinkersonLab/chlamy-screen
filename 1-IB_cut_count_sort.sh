#!/bin/bash

# 1- cutadapt to cut out sequence
# 2- count barcodes
# 3- reverse barcode seq, prep columns

#data location
loc="/shared/Labs/Jonikas/Everyone/deepseq_data/160314_Rnd2_lib_pooled_screens/"

#make dict to pull names of samples (convert names of reads to sample names)
#1st set
declare -A samples=(['160219_HAVERS_0552_BC7VJJACXX_L4_ATCACG_pf.fastq.gz']='TAP-Light_I' ['160219_HAVERS_0552_BC7VJJACXX_L4_TTAGGC_pf.fastq.gz']='TAP-Light_II' ['160219_HAVERS_0552_BC7VJJACXX_L4_TGACCA_pf.fastq.gz']='TAP-Dark_I' ['160224_TENNISON_0408_AC7Y3MACXX_L1_ATCACG_pf.fastq.gz']='TAP-Dark_II' ['160224_TENNISON_0408_AC7Y3MACXX_L1_TTAGGC_pf.fastq.gz']='TAP-Dark_2l' ['160224_TENNISON_0408_AC7Y3MACXX_L1_TGACCA_pf.fastq.gz']='Zeocin' ['160224_TENNISON_0408_AC7Y3MACXX_L2_ATCACG_pf.fastq.gz']='MMS' ['160224_TENNISON_0408_AC7Y3MACXX_L2_TTAGGC_pf.fastq.gz']='MMC' ['160224_TENNISON_0408_AC7Y3MACXX_L2_TGACCA_pf.fastq.gz']='Cis-Platin' ['160224_TENNISON_0408_AC7Y3MACXX_L3_ATCACG_pf.fastq.gz']='TP-light_Air_I' ['160224_TENNISON_0408_AC7Y3MACXX_L3_TTAGGC_pf.fastq.gz']='TP-light_Air_I_2nd' ['160224_TENNISON_0408_AC7Y3MACXX_L3_TGACCA_pf.fastq.gz']='TP-light_Air_II' ['160224_TENNISON_0408_AC7Y3MACXX_L4_ATCACG_pf.fastq.gz']='TP-light_CO2_I' ['160224_TENNISON_0408_AC7Y3MACXX_L4_TTAGGC_pf.fastq.gz']='TP-light_CO2_II' ['160224_TENNISON_0408_AC7Y3MACXX_L4_TGACCA_pf.fastq.gz']='Anoxia' ['160224_TENNISON_0408_AC7Y3MACXX_L6_ATCACG_pf.fastq.gz']='Flagella_D1H_S' ['160224_TENNISON_0408_AC7Y3MACXX_L6_TTAGGC_pf.fastq.gz']='NaCl' ['160224_TENNISON_0408_AC7Y3MACXX_L6_TGACCA_pf.fastq.gz']='Mannitol' ['160224_TENNISON_0408_AC7Y3MACXX_L7_ATCACG_pf.fastq.gz']='N7' ['160224_TENNISON_0408_AC7Y3MACXX_L7_TTAGGC_pf.fastq.gz']='N8' ['160224_TENNISON_0408_AC7Y3MACXX_L7_TGACCA_pf.fastq.gz']='N9' ['160224_TENNISON_0408_AC7Y3MACXX_L8_ATCACG_pf.fastq.gz']='N11' ['160224_TENNISON_0408_AC7Y3MACXX_L8_TTAGGC_pf.fastq.gz']='N12' ['160224_TENNISON_0408_AC7Y3MACXX_L8_TGACCA_pf.fastq.gz']='N13')

for index in ${loc}*.gz

do
echo ""
echo "----------------------------------------------------------------"
echo $index

# 1- cutadapt to cut out sequence
echo R2-${samples["${index#$loc}"]}
## 5 prime IB
#cutadapt -a GGCAAGCTAGAGA -e 0.1 -m 21 -M 23 150807_BRISCOE_0249_BC7CMJACXX_L4_${index}_pf.fastq.gz -o ${index} # 5 prime IB
## 3 prime IB
cutadapt -a TAGCGCGGGGCGT -e 0.1 -m 21 -M 23 ${index} -o R2-${samples["${index#$loc}"]}.fastq # 3 prime IB

# 2- count barcodes
echo R2-${samples["${index#$loc}"]}.count
fastx_collapser -i R2-${samples["${index#$loc}"]}.fastq -o R2-${samples["${index#$loc}"]}.count

# 3- reverse barcode seq, prep columns
echo $index.count.sort
revseq R2-${samples["${index#$loc}"]}.count stdout |awk -F "-" '{getline a; print a, $2}' |sort > R2-${samples["${index#$loc}"]}.count.sort

done
