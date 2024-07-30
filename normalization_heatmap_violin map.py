# -*- coding: utf-8 -*-
"""
Normalizing the seqeunce counts by control sequence data

"""
import pandas as pd
from Bio import motifs
import matplotlib.pyplot as plt
import seaborn as sns

'input: the control and un-normalized sequence counts files'

input_C = 'AU_control.csv'
input_L = 'AU_filtered.csv'

'output: normalized sequence counts; heatmap; violinplot'

output_1 = 'AU_normed.csv'
output_2 = 'AU_normed_heatmap.png'
output_3 = 'AU_normed_violinplot.png'


'Load the input files'
df_C = pd.read_csv(input_C, index_col = 0)
df_L = pd.read_csv(input_L, index_col = 0)


'Calculate the normalization coeffecient, alpha'
total_c = df_C['Counts'].sum()
df_C['alpha'] =(df_C['Counts']/total_c)/(1/(4**4))
df_C.set_index(title, inplace = True)
df_C_alpha = df_C[['alpha']]

'normalizing the loop closing sequencing counts: divide the ligation counts by alpha and save the file'

df_L.set_index(title, inplace = True)
df_L_tonorm = pd.concat([df_L, df_C_alpha], axis=1, join = 'inner')
df_L_tonorm['normalized counts'] = round(df_L_tonorm['Counts']/df_L_tonorm['alpha'])
df_L_tonorm.reset_index(inplace =True)
df_L_tonorm['rank before norm']=df_L_tonorm['Counts'].rank(method='dense', ascending = False)
df_L_tonorm['rank after norm']=df_L_tonorm['normalized counts'].rank(method='dense', ascending = False)
df_L_tonorm['rank diff']= df_L_tonorm['rank after norm']-df_L_tonorm['rank before norm']
df_L_tonorm['normalized percentage'] = df_L_tonorm['normalized counts']/df_L_tonorm['normalized counts'].sum()*100
df_L_tonorm['normalized percentage']=df_L_tonorm['normalized percentage'].round(2)
df_L_tonorm = df_L_tonorm.sort_values(by = 'normalized counts', ascending = False)
df_L_tonorm.to_csv(output_1)

'plot the heatmap and save the file'

df_L_tonorm = df_L_tonorm.replace('U', 'T', regex = True)
Seq_all =[]
    
for i in range(len(df_L_tonorm.index)):
    Seq_extend = [df_L_tonorm[title][i]]*int(df_L_tonorm['normalized counts'][i])
    Seq_all = Seq_all + Seq_extend
    
m = motifs.create(Seq_all)
pwm = m.counts.normalize()
df_pwm = pd.DataFrame(pwm)
df_pwm.index += 1
df_pwm = df_pwm.transpose()
df_pwm = df_pwm.rename(index = {'T':'U'})
     
countplot = sns.heatmap(df_pwm, annot=True, fmt = '.2f', cmap="YlGnBu")
countplot.set_title(title+'-normalized')
countplot.set_xlabel('position of nucleobase')
fig = countplot.get_figure()
fig.savefig(output_2)

plt.clf()

'plot the violinplot'

vplot=sns.violinplot(data=df_L_tonorm[['normalized percentage']])
vplot.set_title(title+'_normalized')
fig_v = vplot.get_figure()
fig_v.savefig(output_3)