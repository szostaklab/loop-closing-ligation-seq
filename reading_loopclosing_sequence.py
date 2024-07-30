# -*- coding: utf-8 -*-
"""
Preprocessing the sequence and get the base counts

This piece of scipt is for reading the loop closing sequence with the close base pairing AU;
For other close base pairing, just need to change the fix region sequence;

"""
import Bio
from Bio import SeqIO
import gzip
import numpy as np
import pandas as pd
import math
import timeit

'input: forward read file R1, reverse read file R2'
filenames_r1 = 'R1.fastq.gz'
filenames_r2 = 'R2.fastq.gz'

'output_1: AU_filtered: csv file summarizing the base counts;'
'output_2: AU_filtered_perct:csv file summarizing the sequence number\
 and percentage before and after data preprocessing'
output_1 = 'AU_filtered.csv'
output_2 = 'AU_filtered_perct.csv'


'function read_seq: "load the r1 and r2 fastq.gz file, read the sequences and q values"'
def read_seq(f1, f2):
        
    record_list_r1 = []
    record_list_r2 = []
    
    record_qlist_r1 = []
    record_qlist_r2 = []
    
    with gzip.open(f1, "rt") as handle:
        for seq_record in SeqIO.parse(handle, "fastq"):
            record_list_r1.append(seq_record.seq)
            record_qlist_r1.append(seq_record.letter_annotations["phred_quality"])
            
    with gzip.open(f2, "rt") as handle:
        for seq_record in SeqIO.parse(handle, "fastq"):
            record_list_r2.append(seq_record.seq)
            record_qlist_r2.append(seq_record.letter_annotations["phred_quality"])
    
    return record_list_r1, record_list_r2, record_qlist_r1, record_qlist_r2

"function seq_filter_1: preprocessing step 1;\
 mean of q in R1>=30 & mean of q in R2 >=30. \
 mean of q = -10log(p_mean), take the intersection part" 
     
def seq_filter_1(record_list_r1, record_list_r2, record_qlist_r1, record_qlist_r2):
          
    f1_index_r1 = []
    f1_index_r2 = []
           
    for index, i in enumerate(record_qlist_r1):
        p=[]
        for k in range(len(i)):
            p.append(10**(-i[k]/10))
        if -10*math.log10(np.mean(p))>=30:
            f1_index_r1.append(index)
          
    for index, j in enumerate(record_qlist_r2): 
        p_r2=[]
        for k in range(len(j)):
            p_r2.append(10**(-j[k]/10))
        if -10*math.log10(np.mean(p_r2))>=30:
            f1_index_r2.append(index)   
             
    index_r1 = set(f1_index_r1)
    intersection = index_r1.intersection(f1_index_r2)
    f1_index_r1_r2 = list(intersection)        
    
    list_f1_r1 = [record_list_r1[index] for index in f1_index_r1_r2]
    list_f1_r2 = [record_list_r2[index] for index in f1_index_r1_r2]
      
    return list_f1_r1, list_f1_r2, f1_index_r1_r2, f1_index_r1, f1_index_r2
     

'function seq_filter_2: prepocessing step 2;\
R1[fix seq] perfect match & R2[fix seq] perfect match& R1[N] complementatry to R2[N]\
the fix region [fix seq] depends on the close base-pairing base' 

def seq_filter_2(list_f1_r1, list_f1_r2, f1_index_r1_r2, N):
          
    f2_index_r1=[]
    f2_index_r2=[]
    list_f3_N=[]
    
    for index, i in enumerate(list_f1_r1):
        if i[0:16] == 'TATACGTAAGCAGCGA' and i[16+N:32+N] == 'TCGCTGCTTACGATAA':
            f2_index_r1.append(index)
        
    for index, j in enumerate(list_f1_r2): 
        if j[0:16] == 'TTATCGTAAGCAGCGA' and j[16+N:32+N] == 'TCGCTGCTTACGTATA':
            f2_index_r2.append(index)                                         
   
    f2_list_1 = set(f2_index_r1)
    intersection = f2_list_1.intersection(f2_index_r2)
    f2_index_r1_r2 = list(intersection)
    
    for k in f2_index_r1_r2:
        if list_f1_r1[k][16:16+N].reverse_complement()==list_f1_r2[k][16:16+N]:
            list_f3_N.append(list_f1_r1[k][16:16+N].transcribe())
        
    return f2_index_r1_r2, list_f3_N
  
N = 4
    
"read the record seq list"
record_list_r1, record_list_r2, record_qlist_r1, record_qlist_r2 =\
    read_seq(filenames_r1, filenames_r2)
    
"Filter 1: mean of q = -10 log(mean(p) in R1 and R2 >=30"
list_f1_r1, list_f1_r2, f1_index_r1_r2, f1_index_r1, f1_index_r2 =\
    seq_filter_1(record_list_r1, record_list_r2, record_qlist_r1, record_qlist_r2)

 
"filter 2: R1[fix seq] perfect match & R2[fix seq] perfect match& R1[N] complementatry to R2[N]"
f2_index_r1_r2, list_f3_N = seq_filter_2(list_f1_r1, list_f1_r2, f1_index_r1_r2, N)


"Generate the output and save the files"
df_N = pd.DataFrame(list_f3_N, columns = ['Seq_N'+str(N)])
count_N = df_N.groupby(['Seq_N'+str(N)]).size().sort_values(ascending=False)
df_N_count = pd.DataFrame(count_N, columns = ['Counts'] ).reset_index()
df_N_count['Percentage'] = df_N_count['Counts']/len(list_f3_N)*100
df_N_count['Percentage'] = df_N_count['Percentage'].round(2)
    
    
df_N_count.to_csv(output_1)


data_seq_number_sum = {'filename': filenames_r1, 'total_seq_number': len(record_list_r1),\
                       'seq_filter1_number':len(list_f1_r1), 'seq_filter2_number': \
                           len(f2_index_r1_r2), 'seq_filter3_number': len(list_f3_N) }

df_seq_number_sum = pd.DataFrame.from_dict(data_seq_number_sum, orient = 'index').transpose()
df_seq_number_sum ['filter_seq_percentage'] = df_seq_number_sum['seq_filter3_number']/df_seq_number_sum['total_seq_number']*100
df_seq_number_sum ['filter_seq_percentage'] = df_seq_number_sum['filter_seq_percentage'].astype(float).round(2)
df_seq_number_sum.set_index('filename', inplace = True)

df_seq_number_sum.to_csv(output_2)




