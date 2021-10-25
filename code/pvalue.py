import pandas as pd
import numpy as np
import scipy.stats as stats
from sys import argv
import collections


script,data,label,classes = argv
classes=int(classes)

gene=pd.read_csv('test.csv')

gnumber=gene.shape[0]
#print(gnumber)
targetid=gene.index
gene=gene.reset_index(drop=True)

gene=gene.T
gene=gene.reset_index(drop=True)

label=pd.read_csv(label)
gene['label']=label

group_label=gene.groupby('label')

cellpvalue=pd.DataFrame()
i=0
#deal with the problem pvalue
while i<gnumber :
	count=0
	j=1#j,k means the class number
	while j<=classes :
		k=j+1
		while k<=classes:
			if True:
				u12,pvalue=stats.mannwhitneyu(group_label.get_group(j)[i],group_label.get_group(k)[i])
				cellpvalue.loc[i,count]=pvalue
				count=count+1
			else:
				cellpvalue.loc[i,count]=0.5
			k=k+1
		j=j+1
	i=i+1

#cellpvalue.to_csv('cellpvalue.csv',sep=',',index=False)

pvalue=pd.DataFrame()
pvalue['targetid']=targetid
pvalue['min_value']=cellpvalue.min(axis=1)
pvalue=pvalue.sort_values(by="min_value" , ascending=True)
pvalue.to_csv('pvalue.csv',sep=',',index=False)





